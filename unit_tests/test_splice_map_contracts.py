"""Splice-map file contracts with rbp-maps isolated behind its public interfaces."""

import os
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest

from conftest import workflow_utils
from test_pipeline_contracts import write_bam


@pytest.fixture
def density_api(monkeypatch):
    density = ModuleType("density")
    density.__path__ = []
    read = ModuleType("density.ReadDensity")
    norm = ModuleType("density.normalization_functions")
    maps = ModuleType("density.Map")
    handles = []

    def rd(**kwargs):
        handle = SimpleNamespace(pos=Mock(), neg=Mock(), bam=Mock())
        handles.append(handle)
        return handle

    read.ReadDensity = Mock(side_effect=rd)
    norm.per_region_subtract_and_normalize = Mock()
    norm.bottom_top_values_from_dataframe = Mock(
        side_effect=lambda frame, *args: (frame.min().tolist(), frame.max().tolist())
    )
    maps.SkippedExon = Mock(return_value=Mock())
    density.ReadDensity = read
    density.normalization_functions = norm
    density.Map = maps
    for name, module in [
        ("density", density),
        ("density.ReadDensity", read),
        ("density.normalization_functions", norm),
        ("density.Map", maps),
    ]:
        monkeypatch.setitem(sys.modules, name, module)
    return SimpleNamespace(read=read, norm=norm, maps=maps, handles=handles)


def test_map_output_names_and_freshness(tmp_path):
    lists = workflow_utils.event_lists("KD", tmp_path / "events")
    for p in lists.values():
        p.parent.mkdir(exist_ok=True)
        p.write_text("events")
        os.utime(p, (100, 100))
    outputs = workflow_utils.map_outputs(tmp_path, "S", "KD", "rep1", lists)
    assert len(outputs) == 9 and len(set(outputs)) == 9
    assert outputs[0] == tmp_path / "KD_KD" / "S_rep1.se.png"
    assert not workflow_utils.map_is_finished(tmp_path, "S", "KD", "rep1", lists)
    for p in outputs:
        p.parent.mkdir(exist_ok=True)
        p.touch()
        os.utime(p, (200, 200))
    assert workflow_utils.map_is_finished(tmp_path, "S", "KD", "rep1", lists)
    os.utime(next(iter(lists.values())), (200, 200))
    assert not workflow_utils.map_is_finished(tmp_path, "S", "KD", "rep1", lists)


@pytest.mark.parametrize("full_counts", [False, True])
def test_splice_map_passes_normalization_annotations_and_closes_handles(
    tmp_path, density_api, full_counts
):
    api = density_api
    ip = dict(pos="ip.pos", neg="ip.neg", bam="ip.bam")
    inp = dict(pos="in.pos", neg="in.neg", bam="in.bam")
    if full_counts:
        ip["mapped"], inp["mapped"] = 1000, 2000
    lists = workflow_utils.event_lists("KD", tmp_path)
    png = tmp_path / "new" / "map.png"
    assert workflow_utils.run_splice_map(ip, inp, lists, png) == png
    args, kwargs = api.maps.SkippedExon.call_args
    assert args[2] == str(png)
    assert list(args[4].items()) == [
        (str(lists[s]), "rmats") for s in workflow_utils.SETS
    ]
    assert kwargs == dict(
        exon_offset=50, intron_offset=300, min_density_threshold=0, conf=0.95
    )
    m = api.maps.SkippedExon.return_value
    assert [call[0] for call in m.method_calls] == [
        "create_matrices",
        "normalize_matrix",
        "create_lines",
        "set_background_and_calculate_significance",
        "write_intermediates_to_csv",
        "plot",
    ]
    m.set_background_and_calculate_significance.assert_called_once_with(
        [str(lists[s]) for s in workflow_utils.SETS[:2]],
        str(lists[workflow_utils.SETS[2]]),
        "permutation",
    )
    for handle in api.handles:
        handle.pos.close.assert_called_once()
        handle.neg.close.assert_called_once()
    if full_counts:
        assert [h.bam.mapped for h in api.handles] == [1000, 2000]
    else:
        for handle in api.handles:
            handle.bam.close.assert_called_once()


def test_load_map_combines_background_extrema_and_finds_exon_trace(
    tmp_path, density_api
):
    lists = workflow_utils.event_lists("KD", tmp_path)
    outputs = workflow_utils.map_outputs(tmp_path, "S", "KD", "rep1", lists)
    outputs[0].parent.mkdir()
    for i, p in enumerate(outputs[-2:]):
        pd.DataFrame([[i - 2, i], [i + 3, i + 5]]).to_csv(p, sep="\t")
    c, strand, start, end = workflow_utils.EXONS["exon13"]["grch38"]
    key = f"id\tgene\tsymbol\t{c}\t{strand}\t{start}\t{end}\trest"
    base = str(outputs[0])[:-4]
    for s in workflow_utils.SETS:
        np.savetxt(f"{base}.{lists[s].name}.means.txt", [0.1, 0.2])
        pd.DataFrame([[0.2, 0.4]], index=[key]).to_csv(
            f"{base}.{lists[s].name}.normed_matrix.txt"
        )
    result = workflow_utils.load_splice_map(tmp_path, "S", "KD", "rep1", lists)
    np.testing.assert_array_equal(result["band"][0], [-2, 0])
    np.testing.assert_array_equal(result["band"][1], [4, 6])
    np.testing.assert_array_equal(result["means"][workflow_utils.SETS[0]], [0.1, 0.2])
    trace, group = workflow_utils.exon_trace(result["matrices"], "exon13")
    np.testing.assert_array_equal(trace, [0.2, 0.4])
    assert group == workflow_utils.SETS[0]
    assert workflow_utils.exon_trace(result["matrices"], "exon24nt") == (None, None)


def test_peak_and_significant_positions_ignore_other_regions():
    trace = np.zeros(1400)
    trace[350 + 267] = 5
    trace[700 + 109] = 7
    assert (
        workflow_utils.peak(trace, "3ss") == -33
        and workflow_utils.peak(trace, "5ss") == 60
    )
    line = np.zeros(1400)
    line[550:635] = 1
    line[750:836] = 1
    assert workflow_utils.positions_above_band(line, np.zeros(1400), "3ss") == 85
    assert workflow_utils.positions_above_band(line, np.zeros(1400), "5ss") == 86
    assert workflow_utils.positions_above_band(line, line, "3ss") == 0


@pytest.mark.parametrize("library", [None, {"mapped_F260": 100, "mapped_records": 200}])
def test_signal_tracks_strand_direction_and_full_library_counts(
    tmp_path, monkeypatch, library
):
    bam = write_bam(tmp_path / "a.bam", [129])
    maps = ModuleType("maps")
    maps.plot_map = SimpleNamespace(ensure_density_bigwigs=Mock())
    monkeypatch.setitem(sys.modules, "maps", maps)
    monkeypatch.setattr(workflow_utils, "read2_only", lambda *a: str(bam))
    monkeypatch.setattr(
        workflow_utils,
        "sense_fraction",
        lambda *a: (0.8 if library is None else 0.2, 10),
    )
    scaled = Mock()
    monkeypatch.setattr(workflow_utils, "_bigwigs_scaled", scaled)
    result = workflow_utils.signal_tracks(
        bam, "bg", tmp_path / "signal", library=library
    )
    assert result["sense_fraction"] == (0.8 if library is None else 0.2)
    assert (tmp_path / "signal" / "a.bam.chrom.sizes").read_text() == "chr1\t10000\n"
    if library is None:
        assert maps.plot_map.ensure_density_bigwigs.call_args.kwargs["direction"] == "f"
        scaled.assert_not_called()
    else:
        assert result["mapped"] == 200 and scaled.call_args.args[-2:] == ("r", 100)
        Path(result["pos"]).touch()
        Path(result["neg"]).touch()
        workflow_utils.signal_tracks(bam, "bg", tmp_path / "signal", library=library)
        assert scaled.call_count == 1


@pytest.mark.parametrize("direction,strands", [("f", ["+", "-"]), ("r", ["-", "+"])])
def test_bigwig_command_scales_rpm_and_sorts_each_strand(
    tmp_path, monkeypatch, direction, strands
):
    processes = []

    def popen(args, **kwargs):
        processes.append(args)
        return SimpleNamespace(stdout=Mock(), wait=lambda: 0)

    monkeypatch.setattr(workflow_utils.subprocess, "Popen", popen)
    convert = Mock()
    monkeypatch.setattr(workflow_utils.subprocess, "check_call", convert)
    workflow_utils._bigwigs_scaled(
        "a.bam", tmp_path / "pos.bw", tmp_path / "neg.bw", "sizes", direction, 100
    )
    genomecov = [args for args in processes if args[0] == "bedtools"]
    assert [a[a.index("-strand") + 1] for a in genomecov] == strands
    assert all(a[-1] == "10000.0" for a in genomecov)
    assert convert.call_count == 2 and not list(tmp_path.glob("*.bg"))
