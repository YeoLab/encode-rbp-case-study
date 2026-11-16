"""Cached HTTP responses and small synthetic ENCODE/RBP-ARK files."""

import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest
import requests

from conftest import workflow_utils, cs3_utils, prepare_events


def archive(path, members):
    with tarfile.open(path, "w:gz") as tar:
        for name, data in members.items():
            data = data.encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return path


@pytest.mark.parametrize("module", [workflow_utils, cs3_utils])
def test_fetch_retries_then_atomically_caches(tmp_path, monkeypatch, module):
    retry = SimpleNamespace(status_code=503, headers={"Retry-After": "2"})
    success = Mock(status_code=200)
    success.iter_content.return_value = [b"abc", b"", b"def"]
    get = Mock(side_effect=[retry, success])
    sleep = Mock()
    monkeypatch.setattr(module.requests, "get", get)
    monkeypatch.setattr(module.time, "sleep", sleep)
    out = tmp_path / "nested" / "file.gz"
    assert module.fetch("https://unit.invalid/file", out, params={"x": 1}) == out
    assert out.read_bytes() == b"abcdef" and not out.with_suffix(".gz.part").exists()
    assert get.call_count == 2
    sleep.assert_called_once_with(2)
    assert module.fetch("unused", out) == out and get.call_count == 2


@pytest.mark.parametrize("module", [workflow_utils, cs3_utils])
def test_fetch_errors_and_zero_attempts(tmp_path, monkeypatch, module):
    response = Mock(status_code=500)
    response.raise_for_status.side_effect = requests.HTTPError("server failure")
    monkeypatch.setattr(module.requests, "get", Mock(return_value=response))
    with pytest.raises(requests.HTTPError):
        module.fetch("unit", tmp_path / "file")
    assert not (tmp_path / "file").exists()
    with pytest.raises(RuntimeError, match="could not download"):
        module.fetch("unit", tmp_path / "file", retries=0)


def test_large_download_checks_size_digest_and_marker(tmp_path, monkeypatch):
    out = tmp_path / "file"
    digest = hashlib.sha256(b"data").hexdigest()

    def curl(args):
        out.write_bytes(b"data")

    call = Mock(side_effect=curl)
    monkeypatch.setattr(workflow_utils.subprocess, "check_call", call)
    assert workflow_utils.download_large("unit", out, 4, digest) == out
    assert workflow_utils.checksum(out, chunk=1) == digest
    assert (
        workflow_utils.download_large("unit", out, 4, digest) == out
        and call.call_count == 1
    )
    out.with_name("file.sha256ok").unlink()
    with pytest.raises(AssertionError, match="checksum mismatch"):
        workflow_utils.download_large("unit", out, 4, "bad")
    monkeypatch.setattr(workflow_utils.subprocess, "check_call", lambda args: None)
    with pytest.raises(AssertionError, match="size"):
        workflow_utils.download_large("unit", out, 5)


def test_json_helpers_send_parameters_and_read_cache(tmp_path, monkeypatch):
    response = Mock()
    response.json.return_value = {"value": 1}
    get = Mock(return_value=response)
    monkeypatch.setattr(cs3_utils.requests, "get", get)
    assert cs3_utils.encode_json("/files/A/", field="x") == {"value": 1}
    assert get.call_args.kwargs["params"] == {"field": "x", "format": "json"}
    payload = tmp_path / "payload.json"
    payload.write_text('{"value": 2}')
    fetch = Mock(return_value=payload)
    monkeypatch.setattr(workflow_utils, "fetch", fetch)
    assert workflow_utils.rbpark_json("/api/gene/", tmp_path, fields="a,b") == {
        "value": 2
    }
    assert fetch.call_args.args[1].name == "api_gene_fields-a+b.json"
    assert workflow_utils.encode_json("/files/A/", tmp_path) == {"value": 2}
    assert fetch.call_args.kwargs["params"] == {"format": "json"}


def test_rbp_ark_records_tables_and_duplicate_files(tmp_path, monkeypatch):
    url = "https://host.invalid/a/file.bam"
    kv = [
        {"key": "path", "value": url},
        {"key": "filesize", "value": "12"},
        {"key": "sha256", "value": "abc"},
    ]
    record = {
        "content": {
            "summary_data": {"basic_stats": {"mapped": [{"key": "IP_1", "value": 20}]}},
            "filelist": [kv, kv],
        }
    }
    api = Mock(return_value=record)
    monkeypatch.setattr(workflow_utils, "rbpark_json", api)
    assert workflow_utils.experiment_record("subject", tmp_path) == record
    assert workflow_utils.read_stats("subject", tmp_path).loc["IP_1", "mapped"] == 20
    result = workflow_utils.file_list("subject", tmp_path)
    assert result.to_dict("records") == [
        dict(url=url, collection_path="/a/file.bam", size=12, sha256="abc")
    ]
    payload = tmp_path / "files.json"
    payload.write_text(
        json.dumps(
            {
                "content": {
                    "filelist": {
                        "a": [kv, {"nested": [kv]}],
                        "ignored": [[{"key": "other"}]],
                    }
                }
            }
        )
    )
    monkeypatch.setattr(cs3_utils, "fetch", lambda *a, **k: payload)
    assert cs3_utils.rbpark_files("subject", tmp_path).url.tolist() == [url]
    monkeypatch.setattr(
        workflow_utils,
        "rbpark_json",
        lambda *a, **k: {
            "experiments": [{"subject": "s", "rbp": "QKI", "windows": [{"start": 1}]}]
        },
    )
    assert workflow_utils.gene_search("G", tmp_path).to_dict("records") == [
        dict(subject="s", rbp="QKI", start=1)
    ]


@pytest.mark.parametrize(
    "module,function",
    [(workflow_utils, "experiments_table"), (cs3_utils, "rbpark_experiments")],
)
def test_experiment_catalog_reads_tsv(tmp_path, monkeypatch, module, function):
    p = tmp_path / "experiments.tsv"
    p.write_text("subject\trbp\ns\tQKI\n")
    monkeypatch.setattr(module, "fetch", lambda *a, **k: p)
    assert getattr(module, function)(tmp_path).rbp.tolist() == ["QKI"]


def test_window_tables_filter_chunks_rename_columns_and_half_open_overlap(
    tmp_path, monkeypatch
):
    p = tmp_path / "windows.tsv.gz"
    table = pd.DataFrame(
        [
            dict(
                subject="s",
                rbp="QKI",
                chr="chr1",
                strand="+",
                gene_name="G",
                gene_id="id",
                start=10,
                end=20,
                enrichment_l2or_mean=3,
                q_max=0.01,
            ),
            dict(
                subject="s",
                rbp="OTHER",
                chr="chr1",
                strand="+",
                gene_name="G",
                gene_id="id",
                start=20,
                end=30,
                enrichment_l2or_mean=2,
                q_max=0.02,
            ),
        ]
    )
    table.to_csv(p, sep="\t", index=False)
    monkeypatch.setattr(cs3_utils, "fetch", lambda *a, **k: p)
    assert len(cs3_utils.rbpark_windows(tmp_path, rbps=["QKI"], chunksize=1)) == 1
    assert len(cs3_utils.rbpark_windows(tmp_path, chunksize=1)) == 2
    monkeypatch.setattr(
        workflow_utils,
        "file_list",
        lambda *a: pd.DataFrame(
            [
                dict(
                    url="unit",
                    collection_path="/reproducible_enriched_windows/s.reproducible_enriched_windows.tsv.gz",
                )
            ]
        ),
    )
    monkeypatch.setattr(workflow_utils, "fetch", lambda *a, **k: p)
    result = workflow_utils.reproducible_windows("s", tmp_path)
    assert result.log2_enrichment.tolist() == [3, 2]
    assert workflow_utils.windows_in(result, "chr1", "+", 10, 20).start.tolist() == [10]
    assert workflow_utils.windows_in(result, "chr1", "-", 10, 20).empty


def test_bam_selection_prefers_read2_and_deduplicates_reused_input(monkeypatch):
    subject = "QKI_K562_5007"
    rows = [
        dict(
            collection_path=f"/bams/dedup/{folder}/{subject}_{sample}.genome.Aligned.sort.dedup{suffix}.bam",
            url="unit",
        )
        for sample in ["IP_1", "IN_1", "IP_2", "IN_2"]
        for folder, suffix in [("genome", ""), ("genome_R2", ".R2")]
    ]
    table = pd.DataFrame(rows)
    assert (
        "genome_R2" in workflow_utils.pick_bam(table, subject, "IP_1").collection_path
    )
    with pytest.raises(AssertionError, match="no deduplicated BAM"):
        workflow_utils.pick_bam(table, subject, "IP_9")
    monkeypatch.setattr(workflow_utils, "MAP_PAIRS", [(subject, "kd")])
    monkeypatch.setattr(workflow_utils, "file_list", lambda *a: table)
    result = workflow_utils.eclip_bams("cache")
    assert sorted(result["sample"]) == ["IN_1", "IN_2", "IP_1", "IP_2"]


def test_encode_knockdown_parents_accept_strings_or_objects(monkeypatch):
    kd = workflow_utils.KNOCKDOWNS["QKI_K562"]

    def metadata(path, cache):
        if path == f"/files/{kd['rmats']}/":
            return {"derived_from": ["/files/a/", {"@id": "/files/b/"}]}
        control = path == "/files/b/"
        return dict(
            dataset={
                "@id": "/experiments/"
                + (kd["control_experiment"] if control else kd["kd_experiment"])
                + "/"
            },
            biological_replicates=[2 if control else 1],
            accession="b" if control else "a",
            assembly="hg19",
            file_size=12,
            md5sum="md5",
            href=path + "bam",
        )

    monkeypatch.setattr(workflow_utils, "encode_json", metadata)
    result = workflow_utils.knockdown_bams("QKI_K562", "cache")
    assert result.condition.tolist() == ["knockdown", "control"]
    assert result.url.str.startswith(workflow_utils.ENCODE).all()
    get = Mock(return_value=Path("out"))
    monkeypatch.setattr(workflow_utils, "download_large", get)
    assert workflow_utils.encode_file("a", "folder", "cache") == Path("out")
    assert get.call_args.args[-1] == "md5"


def test_legacy_knockdown_catalog_and_smallest_successful_tar(monkeypatch, events):
    monkeypatch.setattr(
        cs3_utils,
        "encode_json",
        lambda *a, **k: {
            "@graph": [
                dict(
                    accession="A",
                    dataset="/experiments/E/",
                    biosample_ontology={"term_name": "K562"},
                )
            ]
        },
    )
    files = cs3_utils.knockdown_files("splicing")
    assert files.iloc[0].to_dict() == dict(
        accession="A", experiment="E", target=None, cell="K562", size=0
    )
    files = pd.DataFrame(
        [
            dict(accession="big", target="QKI", cell="K562", size=2),
            dict(accession="small", target="QKI", cell="K562", size=1),
        ]
    )
    calls = []

    def load(acc, cache):
        calls.append(acc)
        return (
            None
            if acc == "small"
            else pd.DataFrame([events(), events(GeneID="ENSG00000151923.2")])
        )

    monkeypatch.setattr(cs3_utils, "rmats_se", load)
    acc, table = cs3_utils.se_for("QKI", "K562", files, "cache")
    assert calls == ["small", "big"] and acc == "big" and len(table) == 1
    assert cs3_utils.se_for("other", "K562", files, "cache") == (None, None)


def test_tar_readers_choose_mats_not_normalized_output(tmp_path, monkeypatch):
    data = 'geneSymbol\tvalue\n"PICALM"\t7\n'
    tar = archive(
        tmp_path / "events.tar.gz",
        {
            "x/MATS_Norm_output/SE.MATS.JunctionCountOnly.txt": "wrong",
            "x/MATS_output/SE.MATS.JunctionCountOnly.txt": data,
        },
    )
    assert prepare_events.se_table(tar) == data.splitlines()
    monkeypatch.setattr(cs3_utils, "fetch", lambda *a, **k: tar)
    assert cs3_utils.rmats_se("A", tmp_path).geneSymbol.tolist() == ["PICALM"]
    out = tmp_path / "extracted" / "se.txt"
    assert workflow_utils.encode_rmats_se(tar, out).read_text() == data
    assert workflow_utils.encode_rmats_se("not-used", out) == out
    bad = archive(tmp_path / "empty.tar.gz", {"other.txt": "empty"})
    with pytest.raises(SystemExit, match="no MATS_output"):
        prepare_events.se_table(bad)
    monkeypatch.setattr(cs3_utils, "fetch", lambda *a, **k: bad)
    assert cs3_utils.rmats_se("A", tmp_path) is None


def test_deseq_tables_match_gene_versions_and_skip_broken_files(tmp_path, monkeypatch):
    p = tmp_path / "table.tsv"
    p.write_text("gene\tlog2FoldChange\nENSG1.2\t-2.5\n")
    assert workflow_utils.deseq2_log2fc(p, "ENSG1") == -2.5
    assert np.isnan(workflow_utils.deseq2_log2fc(p, "missing"))
    monkeypatch.setattr(
        cs3_utils,
        "encode_json",
        lambda *a, **k: {
            "@graph": [
                dict(accession="bad", file_size=1),
                dict(accession="ok", file_size=2),
            ]
        },
    )

    def fetch(url, *a, **k):
        if "/bad/" in url:
            raise requests.HTTPError("not available")
        return p

    monkeypatch.setattr(cs3_utils, "fetch", fetch)
    assert cs3_utils.deseq2_log2fc("E", "ENSG1", tmp_path) == {"ok": -2.5}


@pytest.mark.parametrize("module", [workflow_utils, cs3_utils])
def test_gencode_transcripts_use_zero_based_sorted_exons(tmp_path, monkeypatch, module):
    p = tmp_path / "gencode.gz"
    attrs = 'gene_name "PICALM"; transcript_id "T1"; transcript_name "PICALM-1"; transcript_type "protein_coding"; tag "MANE_Select";'
    with gzip.open(p, "wt") as f:
        f.write("# other gene\n")
        for feature, start, end in [
            ("gene", 1, 99),
            ("transcript", 1, 99),
            ("exon", 31, 40),
            ("exon", 11, 20),
        ]:
            f.write(f"chr1\tG\t{feature}\t{start}\t{end}\t.\t-\t.\t{attrs}\n")
    monkeypatch.setattr(module, "fetch", lambda *a, **k: p)
    tx = module.gencode_transcripts("PICALM", tmp_path)
    assert tx["T1"]["exons"] == [[10, 20], [30, 40]]
    assert tx["T1"]["tags"] == ["MANE_Select"]


@pytest.mark.parametrize("module", [workflow_utils, cs3_utils])
def test_lifter_is_initialized_once_from_cached_chain(tmp_path, monkeypatch, module):
    import pyliftover

    loader = Mock(return_value=object())
    monkeypatch.setattr(pyliftover, "LiftOver", loader)
    monkeypatch.setattr(module, "_LIFT", None)
    monkeypatch.setattr(module, "fetch", lambda *a, **k: tmp_path / "chain.gz")
    assert module.lifter(tmp_path) is module.lifter(tmp_path)
    loader.assert_called_once_with(str(tmp_path / "chain.gz"))


@pytest.mark.parametrize("strand,expected", [("+", "ACGUN"), ("-", "NACGU")])
def test_rna_sequence_both_apis(tmp_path, monkeypatch, strand, expected):
    p = tmp_path / "dna.json"
    p.write_text('{"dna":"acgtn"}')
    monkeypatch.setattr(workflow_utils, "fetch", lambda *a, **k: p)
    response = Mock()
    response.json.return_value = {"dna": "acgtn"}
    monkeypatch.setattr(cs3_utils.requests, "get", Mock(return_value=response))
    assert workflow_utils.rna_sequence("chr1", 0, 5, strand, tmp_path) == expected
    assert cs3_utils.ucsc_sequence("chr1", 0, 5, strand) == expected


@pytest.mark.parametrize("module", [workflow_utils, cs3_utils])
def test_folders_honor_explicit_workspace(tmp_path, monkeypatch, module):
    monkeypatch.setenv("CASE_STUDY_ROOT", str(tmp_path / "root"))
    monkeypatch.setenv("CASE_STUDY_DATA", str(tmp_path / "large"))
    root, data = module.folders()
    assert root == tmp_path / "root" and data == tmp_path / "large"
    assert data.is_dir() and (root / "results").is_dir()
