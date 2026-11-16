"""Offline orchestration, deterministic event lists, and BAM command contracts."""

import gzip
from pathlib import Path
from unittest.mock import Mock

import pandas as pd
import pysam
import pytest

from conftest import workflow_utils, prepare_events
from test_download_contracts import archive
import make_lite_data


def test_event_lists_deduplicate_lift_and_preserve_unchanged_mtime(
    tmp_path, monkeypatch, events, lift
):
    table = pd.DataFrame(
        [
            events(ID="a"),
            events(ID="b", IJC_SAMPLE_1="50,50"),
            events(ID="c", chr="chr2", IncLevelDifference=-0.2),
            events(ID="bg", chr="chr3", FDR=0.8, PValue=0.8),
        ]
    )
    src = tmp_path / "se.tsv"
    table.to_csv(src, sep="\t", index=False)
    monkeypatch.setattr(workflow_utils, "lifter", lambda cache: lift)
    counts, paths = workflow_utils.build_event_lists(
        "KD", src, tmp_path / "lists", tmp_path
    )
    assert counts[workflow_utils.SETS[0]] == {
        "rMATS events": 2,
        "one per exon": 1,
        "lifted to GRCh38": 1,
    }
    result = pd.read_csv(paths[workflow_utils.SETS[0]], sep="\t")
    assert result.ID.tolist() == ["b"] and result.exonStart_0base.tolist() == [1100]
    times = {s: p.stat().st_mtime_ns for s, p in paths.items()}
    assert (
        workflow_utils.build_event_lists("KD", src, tmp_path / "lists", tmp_path)[1]
        == paths
    )
    assert {s: p.stat().st_mtime_ns for s, p in paths.items()} == times
    assert workflow_utils.event_lists("KD", tmp_path / "lists") == paths
    assert workflow_utils.which_set(paths, "exon13") is None
    row = events(**dict(workflow_utils.EXONS["exon13"]["hg19"]))
    c, strand, a, b = workflow_utils.EXONS["exon13"]["grch38"]
    row.update(chr=c, strand=strand, exonStart_0base=a, exonEnd=b)
    pd.DataFrame([row]).to_csv(paths[workflow_utils.SETS[1]], sep="\t", index=False)
    assert workflow_utils.which_set(paths, "exon13") == workflow_utils.SETS[1]


def test_exon_row_exact_and_highest_read_fallback(tmp_path, events):
    info = workflow_utils.EXONS["exon13"]["hg19"]
    exact = events(**info)
    p = tmp_path / "se.tsv"
    pd.DataFrame([exact]).to_csv(p, sep="\t", index=False)
    row = workflow_utils.exon_row(p, "exon13")
    assert row["note"] == "" and workflow_utils.exon_stats(row)["IJC_kd"] == [10, 20]
    assert workflow_utils.exon_stats(row)["psi_kd"] == [0.7, 0.9]
    high = dict(exact, ID="best", upstreamEE=10, IJC_SAMPLE_1="100,200")
    pd.DataFrame([dict(exact, upstreamEE=9), high]).to_csv(p, sep="\t", index=False)
    row = workflow_utils.exon_row(p, "exon13")
    assert row["ID"] == "best" and "flanking exons differ" in row["note"]
    with pytest.raises(AssertionError, match="not in"):
        workflow_utils.exon_row(p, "exon24nt")


def test_isoform_inclusion_and_skipping(locus):
    table = workflow_utils.isoform_table(locus["transcripts"]).set_index("name")
    assert table.loc["PICALM-1", "150-nt exon included"]
    assert table.loc["PICALM-1", "MANE"]
    assert table.loc["PICALM-2", "150-nt exon skipped"]
    assert not table.loc["PICALM-3", "150-nt exon skipped"]


def manifest():
    return pd.DataFrame(
        [
            dict(
                knockdown="KD",
                rep=rep,
                condition=condition,
                bam=f"{condition}{rep}.bam",
            )
            for rep in [2, 1]
            for condition in ["control", "knockdown"]
        ]
    )


def test_rmats_command_and_reuse_invalidates_changed_settings(tmp_path, monkeypatch):
    calls = []

    def run(cmd, **kwargs):
        calls.append(cmd)
        od = Path(cmd[cmd.index("--od") + 1])
        (od / "SE.MATS.JC.txt").write_text("result")
        return 0

    monkeypatch.setattr(workflow_utils.subprocess, "call", run)
    out = tmp_path / "runs"
    result, status = workflow_utils.run_rmats(
        "KD", manifest(), "genes.gtf", out, 3, ["python", "rmats.py"], "settings"
    )
    assert status == "ran" and result.read_text() == "result"
    assert (out / "KD" / "b1.txt").read_text() == "knockdown1.bam,knockdown2.bam\n"
    assert (out / "KD" / "b2.txt").read_text() == "control1.bam,control2.bam\n"
    cmd = calls[0]
    assert (
        cmd[:2] == ["python", "rmats.py"]
        and cmd[cmd.index("--libType") + 1] == "fr-firststrand"
    )
    assert cmd[cmd.index("--nthread") + 1] == "3"
    assert (
        workflow_utils.run_rmats(
            "KD", manifest(), "genes.gtf", out, 3, ["rmats.py"], "settings"
        )[1]
        == "reused"
    )
    stale = out / "KD" / "old"
    stale.mkdir()
    (stale / "file").write_text("old")
    assert (
        workflow_utils.run_rmats(
            "KD", manifest(), "genes.gtf", out, 3, ["rmats.py"], "changed"
        )[1]
        == "ran"
    )
    assert not stale.exists() and len(calls) == 2
    monkeypatch.setattr(workflow_utils.subprocess, "call", lambda *a, **k: 1)
    with pytest.raises(AssertionError, match="rMATS failed"):
        workflow_utils.run_rmats(
            "KD", manifest(), "genes.gtf", out, 3, ["rmats.py"], "failure"
        )
    with pytest.raises(AssertionError):
        workflow_utils.run_rmats("KD", manifest().iloc[:2], "g", out, 1, ["r"], "bad")


def write_bam(path, flags):
    with pysam.AlignmentFile(
        path,
        "wb",
        header={
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 10000}],
        },
    ) as bam:
        for i, flag in enumerate(flags):
            read = pysam.AlignedSegment()
            read.query_name = f"r{i}"
            read.query_sequence = "A" * 20
            read.flag = flag
            read.reference_id = 0
            read.reference_start = 105 + i
            read.mapping_quality = 60
            read.cigar = [(0, 20)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 20)
            bam.write(read)
    pysam.index(str(path))
    return path


def test_bam_read_detection_and_sense_fraction(tmp_path, events):
    bam = write_bam(tmp_path / "paired.bam", [65, 129, 145, 256, 2048])
    assert workflow_utils.has_read1(bam)
    r2 = write_bam(tmp_path / "r2.bam", [129, 145])
    assert not workflow_utils.has_read1(r2, n=0) and not workflow_utils.has_read1(r2)
    p = tmp_path / "events.tsv"
    pd.DataFrame([events(), events(chr="missing")]).to_csv(p, sep="\t", index=False)
    fraction, total = workflow_utils.sense_fraction(bam, p)
    assert total == 3 and fraction == pytest.approx(2 / 3)
    empty = write_bam(tmp_path / "empty.bam", [])
    assert workflow_utils.sense_fraction(empty, p) == (0, 0)


def test_read2_copy_is_created_once(tmp_path, monkeypatch):
    original = tmp_path / "a.bam"
    monkeypatch.setattr(workflow_utils, "has_read1", lambda p: False)
    assert workflow_utils.read2_only(original) == str(original)
    monkeypatch.setattr(workflow_utils, "has_read1", lambda p: True)

    def run(args):
        if args[1] == "view":
            Path(args[args.index("-o") + 1]).touch()

    call = Mock(side_effect=run)
    monkeypatch.setattr(workflow_utils.subprocess, "check_call", call)
    out = workflow_utils.read2_only(original, threads=2)
    assert out.endswith(".R2only.bam") and call.call_count == 2
    assert workflow_utils.read2_only(original) == out and call.call_count == 2
    assert call.call_args_list[0].args[0][2:7] == ["-@", "2", "-f", "128", "-b"]


def test_small_bam_helpers_keep_read2_and_report_full_counts(tmp_path, monkeypatch):
    monkeypatch.setattr(make_lite_data.shutil, "which", lambda cmd: None)
    source = write_bam(tmp_path / "full.bam", [65, 129, 145])
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t0\t500\n")
    out = tmp_path / "lite.bam"
    counts = make_lite_data.lite_eclip(source, bed, out)
    assert counts == dict(
        paired_input=True, mapped_F260=2, mapped_records=2, kept_reads=2
    )
    assert not workflow_utils.has_read1(out)
    assert len(make_lite_data.sh("samtools", "view", "-c", out, capture=True)) > 0
    assert make_lite_data.md5(out) == workflow_utils.checksum(out, "md5")
    normal = Mock(return_value=" output \n")
    monkeypatch.setattr(make_lite_data.subprocess, "check_output", normal)
    assert make_lite_data.sh("echo", "x", capture=True) == "output"
    call = Mock()
    monkeypatch.setattr(make_lite_data.subprocess, "check_call", call)
    make_lite_data.sh("echo", 1)
    call.assert_called_once_with(["echo", "1"])


def test_full_bam_and_bed_helpers(tmp_path):
    (tmp_path / "S_IP_1.genome.Aligned.sort.dedup.R2only.bam").touch()
    full = tmp_path / "S_IP_1.genome.Aligned.sort.dedup.bam"
    full.touch()
    assert make_lite_data.full_eclip_bam(tmp_path, "S", "IP_1") == full
    with pytest.raises(AssertionError, match="no full BAM"):
        make_lite_data.full_eclip_bam(tmp_path, "missing", "IP_1")
    bed = tmp_path / "out.bed"
    assert (
        make_lite_data.write_bed(
            [("chr2", 0, 5), ("chr1", 20, 30), ("chr1", 0, 20)], bed
        )
        == 2
    )
    assert bed.read_text() == "chr1\t0\t30\nchr2\t0\t5\n"
    assert make_lite_data.write_bed([], bed) == 0 and bed.read_text() == ""


def test_splice_regions_use_every_selected_event_and_single_point_lifts(
    tmp_path, monkeypatch, events, lift
):
    p = tmp_path / "se.tsv"
    pd.DataFrame([events()]).to_csv(p, sep="\t", index=False)
    monkeypatch.setattr(workflow_utils, "lifter", lambda cache: lift)
    rows = make_lite_data.splice_site_regions(p, tmp_path)
    assert ("chr1", 700, 1500) in rows and make_lite_data.PICALM_GRCH38 in rows
    monkeypatch.setattr(lift, "convert_coordinate", lambda *a: [])
    assert make_lite_data.splice_site_regions(p, tmp_path) == [
        make_lite_data.PICALM_GRCH38
    ]


def test_picalm_gtf_and_excluded_gene_copy(tmp_path, events):
    gtf = tmp_path / "genes.gtf"
    text = 'chr1\tG\tgene\t5000\t6000\t.\t-\t.\tgene_name "PICALM";\n'
    gtf.write_text(text + text.replace("PICALM", "OTHER"))
    assert make_lite_data.picalm_hg19(gtf) == (("chr1", 3000, 8000), [text])
    src = tmp_path / "all.tsv"
    dst = tmp_path / "kept.tsv"
    pd.DataFrame([events(), events(GeneID="ENSG00000151923.3")]).to_csv(
        src, sep="\t", index=False
    )
    make_lite_data.copy_without_excluded(src, dst)
    assert len(pd.read_csv(dst, sep="\t")) == 1


def test_download_event_main_writes_merged_regions_and_counts(
    tmp_path, monkeypatch, events
):
    chain = tmp_path / "chain.gz"
    with gzip.open(chain, "wt") as f:
        f.write("chain 1 chr1 10000 + 0 10000 chr1 10000 + 0 10000 1\n10000\n")
    table = pd.DataFrame(
        [
            events(),
            events(ID="bg", FDR=0.8),
            events(ID="skip", GeneID="ENSG00000151923"),
        ]
    )
    archive(
        tmp_path / "A.tar.gz",
        {
            "x/MATS_output/SE.MATS.JunctionCountOnly.txt": table.to_csv(
                sep="\t", index=False
            )
        },
    )
    monkeypatch.setattr(prepare_events, "KNOCKDOWNS", {"KD": "A"})
    monkeypatch.setenv("CS3_SKIP_SELFTEST", "1")
    out = tmp_path / "out"
    prepare_events.main(tmp_path, chain, out)
    assert pd.read_csv(
        out / "KD.GRCh38.included-upon-knockdown.txt", sep="\t"
    ).ID.tolist() == ["first"]
    assert (out / "regions.bed").read_text().startswith("chr1\t0\t600\n")
    assert "KD\tincluded-upon-knockdown\t1\t1" in (out / "event_counts.tsv").read_text()
    monkeypatch.delenv("CS3_SKIP_SELFTEST")
    with pytest.raises(SystemExit, match="liftover self-test failed"):
        prepare_events.main(tmp_path, chain, out)


def test_lite_paths_and_library_sizes(tmp_path, monkeypatch):
    monkeypatch.setattr(workflow_utils, "LITE", False)
    assert workflow_utils.data_dirs(tmp_path)["rmats"] == tmp_path / "rmats" / "rerun"
    monkeypatch.setattr(workflow_utils, "LITE", True)
    monkeypatch.setattr(workflow_utils, "LITE_DIR", tmp_path)
    assert workflow_utils.data_dirs(tmp_path)["rmats"] == tmp_path / "rmats_full"
    (tmp_path / "eclip").mkdir()
    (tmp_path / "eclip" / "library_sizes.tsv").write_text(
        "file\tmapped_F260\na.bam\t10\n"
    )
    assert workflow_utils.lite_eclip_bams().bam.tolist() == [
        str(tmp_path / "eclip" / "a.bam")
    ]
    assert workflow_utils.lite_shrna_bam("A") == tmp_path / "shrna" / "A.PICALM.bam"
