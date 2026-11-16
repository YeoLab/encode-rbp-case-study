"""Notebook creation, synthetic figure rendering, and installer isolation."""

import runpy
import shutil
import sys
from pathlib import Path
from unittest.mock import Mock

import matplotlib as mpl
import matplotlib.pyplot as plt
import nbformat
import numpy as np
import pytest

from conftest import ROOT, workflow_utils
import figure
import build_notebooks
import colab_setup


def test_workflow_figure_writes_formats_and_restores_style(
    tmp_path, locus, knockdown, monkeypatch
):
    line = np.linspace(0.0001, 0.0002, 1400)
    maps = [
        dict(
            label=label,
            means={
                s: line * scale for s, scale in zip(workflow_utils.SETS, [1.2, 0.8, 1])
            },
            band=(line * 0.9, line * 1.1),
            traces={"exon13": line, "exon24nt": line[::-1]},
        )
        for label in ["K562", "HepG2"]
    ]
    before = dict(mpl.rcParams)
    saved = []
    close = plt.close
    monkeypatch.setattr(
        plt, "close", lambda f: saved.append(f) if hasattr(f, "axes") else close(f)
    )
    out = tmp_path / "workflow"
    try:
        assert (
            figure.draw_figure(locus, knockdown, maps, out, "synthetic fixtures")
            == str(out) + ".png"
        )
        assert dict(mpl.rcParams) == before
        assert out.with_suffix(".png").read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
        assert out.with_suffix(".pdf").read_bytes().startswith(b"%PDF")
        texts = [text.get_text() for ax in saved[0].axes for text in ax.texts] + [
            text.get_text() for text in saved[0].texts
        ]
        assert set("abcde").issubset(texts)
        assert any("synthetic fixtures" in text for text in texts)
        assert any("PICALM-1" in text for text in texts)
        assert {str(number) for number in range(1, 8)}.issubset(texts)
    finally:
        for f in saved:
            close(f)


def assert_clipping_warning(warn, capsys):
    fig = plt.figure(figsize=(3, 3))
    try:
        fig.add_axes([0.15, 0.2, 0.5, 0.5], title="inside axes")
        fig.text(0.5, 0.9, "inside label", ha="center")
        warn(fig)
        assert capsys.readouterr().out == ""

        fig.add_axes([0.8, 0.2, 0.4, 0.5], title="clipped axes")
        fig.text(-0.25, 0.95, "clipped label")
        warn(fig)
        assert capsys.readouterr().out.splitlines() == [
            "warning: cut off at the figure edge: 'clipped axes'",
            "warning: cut off at the figure edge: 'clipped label'",
        ]
    finally:
        plt.close(fig)


def test_workflow_figure_warns_when_artists_extend_beyond_canvas(capsys):
    assert_clipping_warning(figure.warn_if_cut_off, capsys)


def test_standalone_redraw_uses_bundled_small_inputs(tmp_path, monkeypatch, capsys):
    src = ROOT / "data_and_code" / "picalm"
    for name in ["picalm.json", "picalm_locus.json"]:
        shutil.copyfile(src / name, tmp_path / name)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("MAP_DIR", str(src / "maps_stranded"))
    monkeypatch.setenv("MAP_NOTE", "offline redraw")
    with plt.rc_context():
        try:
            ns = runpy.run_path(str(src / "make_fig.py"))
            assert (
                (tmp_path / "CS3_PICALM_QKI_draft.pdf").read_bytes().startswith(b"%PDF")
            )
            assert (tmp_path / "CS3_PICALM_QKI_draft.png").stat().st_size > 10000
            assert ns["picalm_trace"](
                "QKI_K562_5007", "QKI_HepG2", "rep1", "exon13"
            ).shape == (1400,)
            assert any("offline redraw" in text.get_text() for text in ns["fig"].texts)
            capsys.readouterr()
            assert_clipping_warning(ns["warn_if_cut_off"], capsys)
        finally:
            plt.close("all")


def test_psi_chart_groups_conditions(knockdown):
    fig = workflow_utils.plot_psi(knockdown)
    try:
        assert len(fig.axes) == 2
        assert all(len(ax.patches) == 4 for ax in fig.axes)
        assert [p.get_height() for p in fig.axes[0].patches] == pytest.approx(
            [0.35, 0.65, 0.35, 0.65]
        )
    finally:
        plt.close(fig)


def test_notebook_builder_embeds_screenshots_and_compiles_cells(tmp_path):
    build_notebooks.build(tmp_path)
    assert {p.name for p in tmp_path.glob("*.ipynb")} == set(
        build_notebooks.TITLES.values()
    )
    attachments = 0
    for key, name in build_notebooks.TITLES.items():
        nb = nbformat.read(tmp_path / name, as_version=4)
        nbformat.validate(nb)
        assert nb.metadata.kernelspec.name == "python3"
        setup = nb.cells[2].source
        assert "colab_setup.setup" in setup
        assert (", rmats=True" in setup) == (
            key in ["05_differential_splicing", "00_colab_run_all"]
        )
        for i, cell in enumerate(nb.cells):
            if cell.cell_type == "code":
                compile(cell.source, f"{key}:{i}", "exec")
            for attachment in cell.get("attachments", {}).values():
                attachments += 1
                assert attachment.get("image/jpeg") or attachment.get("image/png")
    assert attachments >= 8
    assert build_notebooks.shot("s.png", "description") == (
        "shot",
        "s.png",
        "description",
    )


def test_colab_installers_pin_core_and_reuse_success_marker(tmp_path, monkeypatch):
    monkeypatch.setattr(colab_setup, "TOOLS", tmp_path / "tools")
    shell = Mock()
    monkeypatch.setattr(colab_setup, "_sh", shell)
    versions = colab_setup._core_versions()
    assert set(versions) == set(colab_setup.CORE) and all(versions.values())
    colab_setup.install_tools(dry_run=True)
    assert not (colab_setup.TOOLS / ".done").exists()
    assert (colab_setup.TOOLS / "constraints.txt").read_text() == "".join(
        f"{p}=={v}\n" for p, v in versions.items()
    )
    assert shell.call_count == 4 and all(
        call.args[1] is True for call in shell.call_args_list
    )
    colab_setup.install_tools()
    assert (colab_setup.TOOLS / ".done").read_text() == "ok\n"
    before = shell.call_count
    colab_setup.install_tools()
    assert shell.call_count == before


def test_shell_helper_dry_run_never_executes(monkeypatch):
    command = Mock()
    monkeypatch.setattr(colab_setup.subprocess, "run", command)
    colab_setup._sh("echo offline", dry_run=True)
    command.assert_not_called()
    colab_setup._sh("echo offline")
    command.assert_called_once_with("echo offline", shell=True, check=True)


def test_rmats_installer_sets_command_and_reuses_environment(tmp_path, monkeypatch):
    monkeypatch.setattr(colab_setup, "TOOLS", tmp_path / "tools")
    monkeypatch.setattr(colab_setup, "RMATS_ENV", tmp_path / "rmats")
    shell = Mock()
    monkeypatch.setattr(colab_setup, "_sh", shell)
    monkeypatch.setenv("RMATS_CMD", "old")
    colab_setup.install_rmats(dry_run=True)
    assert shell.call_count == 2 and "rmats=4.3.0" in shell.call_args.args[0]
    assert (
        str(colab_setup.RMATS_ENV / "bin" / "python")
        in colab_setup.os.environ["RMATS_CMD"]
    )
    script = colab_setup.RMATS_ENV / "bin" / "rmats.py"
    script.parent.mkdir(parents=True)
    script.touch()
    colab_setup.install_rmats()
    assert shell.call_count == 2


@pytest.mark.parametrize("lite", [True, False])
def test_colab_setup_redirects_outputs_and_preserves_workflow(
    tmp_path, monkeypatch, lite
):
    workflow = tmp_path / "workflow"
    workflow.mkdir()
    (workflow / "workflow_utils.py").touch()
    (workflow / "lite_data").mkdir()
    (workflow / "lite_data" / "MANIFEST.tsv").touch()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "path", list(sys.path))
    for key in [
        "CS3_LITE",
        "CS3_LITE_DIR",
        "CASE_STUDY_ROOT",
        "CASE_STUDY_DATA",
        "PATH",
    ]:
        monkeypatch.setenv(key, colab_setup.os.environ.get(key, ""))
    tools = Mock()
    rmats = Mock()
    limit = Mock()
    monkeypatch.setattr(colab_setup, "install_tools", tools)
    monkeypatch.setattr(colab_setup, "install_rmats", rmats)
    monkeypatch.setattr(colab_setup.resource, "getrlimit", lambda *a: (256, 20000))
    monkeypatch.setattr(colab_setup.resource, "setrlimit", limit)
    output = tmp_path / "output"
    colab_setup.setup(
        workflow,
        lite=lite,
        data_on_drive=True,
        rmats=True,
        output_folder=output,
        dry_run=True,
    )
    suffix = "_lite" if lite else ""
    assert Path(colab_setup.os.environ["CASE_STUDY_ROOT"]) == output / (
        "colab_output" + suffix
    )
    assert Path(colab_setup.os.environ["CASE_STUDY_DATA"]) == output / (
        "colab_data" + suffix
    )
    assert (colab_setup.os.environ.get("CS3_LITE") == "1") is lite
    limit.assert_called_once_with(colab_setup.resource.RLIMIT_NOFILE, (10240, 20000))
    tools.assert_called_once_with(True)
    rmats.assert_called_once_with(dry_run=True)


def test_colab_setup_rejects_missing_workflow(tmp_path):
    with pytest.raises(AssertionError, match="workflow_utils.py not found"):
        colab_setup.setup(tmp_path, dry_run=True)
