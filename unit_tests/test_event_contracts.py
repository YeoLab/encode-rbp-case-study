"""Coordinate boundaries and the distinct published event-selection rules."""

import gzip
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from conftest import workflow_utils, cs3_utils, prepare_events


def as_fields(row):
    return [str(v) for v in row.values()], {k: i for i, k in enumerate(row)}


@pytest.mark.parametrize("module", [workflow_utils, cs3_utils])
def test_gene_ids_remove_quotes_and_version(module):
    assert module.gene_ids(
        pd.DataFrame({"GeneID": ['"ENSG123.4"', "ENSG456"]})
    ).tolist() == ["ENSG123", "ENSG456"]


@pytest.mark.parametrize("module", [workflow_utils, cs3_utils])
def test_one_event_per_exon_merges_touching_spans_and_keeps_stable_ties(module, events):
    rows = [
        events(ID="first"),
        events(ID="tie", strand="-"),
        events(ID="best", upstreamEE=200, downstreamES=400, IJC_SAMPLE_1="40,40"),
        events(ID="separate", upstreamEE=401, downstreamES=500),
        events(ID="chr2", chr="chr2"),
    ]
    fields, col = as_fields(rows[0])
    assert module.one_event_per_exon([as_fields(r)[0] for r in rows], col) == [
        "best",
        "separate",
        "chr2",
    ]
    assert module.one_event_per_exon([fields, as_fields(rows[1])[0]], col) == ["first"]
    assert module.one_event_per_exon([], col) == []


@pytest.mark.parametrize(
    "changes,expected",
    [
        ({}, 0),
        ({"IncLevelDifference": -0.2}, 1),
        ({"FDR": 0.8, "PValue": 0.8}, 2),
        ({"FDR": 0.1}, None),
        ({"IncLevelDifference": 0.05}, None),
        ({"PValue": 0.05}, None),
        ({"FDR": 0.8, "PValue": 0.5}, None),
        ({"FDR": 0.8, "PValue": 0.8, "IncLevel1": "0.05,0.05"}, None),
        ({"FDR": 0.8, "PValue": 0.8, "IncLevel1": "NA,"}, None),
        ({"GeneID": '"ENSG00000151923.9"'}, None),
    ],
)
def test_workflow_classification_strict_thresholds(events, changes, expected):
    masks = workflow_utils.classify(pd.DataFrame([events(**changes)]))
    assert [s for s in workflow_utils.SETS if masks[s].iloc[0]] == (
        [] if expected is None else [workflow_utils.SETS[expected]]
    )


@pytest.mark.parametrize(
    "changes,expected",
    [
        ({}, 0),
        ({"IncLevelDifference": -0.05}, 1),
        ({"FDR": 0.5}, 2),
        (
            {
                "IJC_SAMPLE_1": "0",
                "IJC_SAMPLE_2": "0",
                "SJC_SAMPLE_1": "0",
                "SJC_SAMPLE_2": "0",
            },
            None,
        ),
        ({"GeneID": '"ENSG00000151923.9"'}, None),
        ({"FDR": 0.05}, None),
        ({"FDR": 0.8, "IncLevel2": "NA,"}, None),
    ],
)
def test_download_event_classification_read_floor(events, changes, expected):
    f, col = as_fields(events(**changes))
    assert prepare_events.classify(f, col) == (
        None if expected is None else prepare_events.SETS[expected]
    )


def test_legacy_event_sets_include_boundary_and_do_not_mutate(events):
    table = pd.DataFrame(
        [
            events(ID="inc", IncLevelDifference=0.05),
            events(ID="exc", IncLevelDifference=-0.05),
            events(ID="bg", FDR=0.5),
            events(ID="missing", IncLevel1="NA,", IncLevel2="NA,"),
        ]
    )
    result = cs3_utils.event_sets(table)
    assert result["inc"].ID.tolist() == ["inc", "missing"]
    assert result["exc"].ID.tolist() == ["exc"]
    assert result["bg"].ID.tolist() == ["bg"]
    stats = cs3_utils.add_se_stats(table)
    assert stats.junction_reads.iloc[0] == 65
    assert stats.psi_kd.iloc[0] == pytest.approx(0.8)
    assert np.isnan(stats.psi_ctl.iloc[-1])
    assert "junction_reads" not in table


@pytest.mark.parametrize(
    "text,expected", [("1,NA,2,", [1, 2]), ("NA,", []), ("2.5", [2.5])]
)
def test_download_numeric_fields(text, expected):
    assert prepare_events.nums(text) == expected


@pytest.mark.parametrize(
    "mode",
    ["valid", "missing", "ambiguous", "reverse", "length", "chromosome", "order"],
)
def test_workflow_liftover_requires_whole_ordered_unique_exons(events, mode):
    f, col = as_fields(events())

    def convert(chrom, pos):
        if mode == "missing" and pos == 100:
            return []
        if mode == "ambiguous" and pos == 100:
            return [(chrom, 1100, "+"), (chrom, 2100, "+")]
        target = "chr2" if mode == "chromosome" and pos < 50 else chrom
        p = pos + (0 if mode == "order" and pos >= 200 else 1000)
        if mode == "length" and pos == 119:
            p += 1
        return [(target, p, "-" if mode == "reverse" else "+")]

    out = workflow_utils.lift_event(f, col, SimpleNamespace(convert_coordinate=convert))
    if mode == "valid":
        assert out[col["exonStart_0base"]] == "1100"
        assert out[col["exonEnd"]] == "1120"
        assert f[col["exonStart_0base"]] == "100"
    else:
        assert out is None


def test_chain_coordinates_gaps_reverse_and_unknown(tmp_path, events):
    path = tmp_path / "chain.gz"
    with gzip.open(path, "wt") as stream:
        stream.write(
            "chain 100 chr1 1000 + 0 50 chr1 2000 + 100 150 1\n20 10 5\n20\n\n"
            "chain 50 chr2 1000 + 0 10 chr3 2000 - 100 110 2\n10\n"
        )
    chain = prepare_events.Chain(path)
    assert chain.convert("chr1", 0) == [("chr1", 100, "+", 100)]
    assert chain.convert("chr1", 19)[0][1] == 119
    assert chain.convert("chr1", 20) == []
    assert chain.convert("chr1", 30)[0][1] == 125
    assert chain.convert("chr2", 0) == [("chr3", 1899, "-", 50)]
    assert chain.convert("missing", 0) == []
    assert prepare_events.lift_point(chain, "chr1", "0") == ("chr1", 100, "+")
    assert prepare_events.lift_point(chain, "chr1", 20) is None


@pytest.mark.parametrize("mode", ["valid", "missing", "length", "reverse", "other"])
def test_download_lift_fields_preserves_length(events, mode):
    f, col = as_fields(events())

    def convert(chrom, pos):
        if mode == "missing":
            return []
        return [
            (
                "chr2" if mode == "other" else chrom,
                pos + 1000 + (1 if mode == "length" and pos == 119 else 0),
                "-" if mode == "reverse" else "+",
                100,
            )
        ]

    out = prepare_events.lift_fields(SimpleNamespace(convert=convert), f, col)
    if mode == "valid":
        assert out[col["exonEnd"]] == "1120"
    else:
        assert out is None


def test_legacy_lift_events_rejects_missing_wrong_chromosome_and_resizing(
    monkeypatch, events, lift
):
    monkeypatch.setattr(cs3_utils, "lifter", lambda cache: lift)
    table = pd.DataFrame([events()])
    assert cs3_utils.lift_events(table, "cache").exonEnd_hg38.tolist() == [1120]
    for hits in [[], [("other", 100, "+", 1)], [("chr1", 100, "+", 1)]]:
        monkeypatch.setattr(lift, "convert_coordinate", lambda *args: hits)
        assert cs3_utils.lift_events(table, "cache").empty


def test_window_index_is_half_open_strand_specific_and_filters_quality():
    windows = pd.DataFrame(
        [
            dict(
                subject="QKI",
                chr="chr1",
                strand="+",
                start=100,
                end=120,
                q_value=0.05,
                log2_enrichment=2,
            ),
            dict(
                subject="QKI",
                chr="chr1",
                strand="+",
                start=200,
                end=220,
                q_value=0.1,
                log2_enrichment=4,
            ),
            dict(
                subject="QKI",
                chr="chr1",
                strand="+",
                start=300,
                end=320,
                q_value=0.01,
                log2_enrichment=1,
            ),
        ]
    )
    idx = cs3_utils.WindowIndex(windows)
    assert idx.any("QKI", "chr1", "+", 119, 121)
    for a, b in [(80, 100), (120, 140), (200, 220), (300, 320)]:
        assert not idx.any("QKI", "chr1", "+", a, b)
    assert not idx.any("QKI", "chr1", "-", 100, 120)
    events = pd.DataFrame(
        [dict(chr="chr1", strand="+", exonStart_0base_hg38=100, exonEnd_hg38=120)]
    )
    assert cs3_utils.window_map(events, idx, "QKI").tolist() == [0] * 6 + [1] + [0] * 6


@pytest.mark.parametrize(
    "strand,first,last",
    [("+", (700, 750), (1270, 1320)), ("-", (1270, 1320), (700, 750))],
)
def test_exon_bins_follow_transcription(strand, first, last):
    bins = cs3_utils.exon_bins(strand, 1000, 1020)
    assert len(bins) == 13 and bins[6] == (1000, 1020)
    assert bins[0] == first and bins[-1] == last


@pytest.mark.parametrize("strand,left,right", [("+", -1, 1), ("-", 1, -1)])
def test_relative_position_counts_intronic_nucleotides(
    monkeypatch, strand, left, right
):
    monkeypatch.setattr(
        workflow_utils, "EXONS", {"test": {"grch38": ("chr1", strand, 100, 120)}}
    )
    assert workflow_utils.relative_position(90, 100, "test") == left
    assert workflow_utils.relative_position(120, 130, "test") == right
    assert workflow_utils.relative_position(99, 101, "test") == 0
