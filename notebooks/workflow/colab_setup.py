"""Set up a Google Colab session for the Case Study 3 workflow notebooks.

Each notebook's first cell calls setup() when it runs on Colab (and does nothing elsewhere). That cell clones the
repository onto the Colab machine (or uses a copy on Google Drive). setup() installs the command-line tools and Python
packages once per session, and points the notebooks at two folders:

    results and the figure    <output folder>/colab_output/ (lite: colab_output_lite/)
    downloads and big files   full: /content/cs3_data/ on the Colab machine, deleted when the session ends
                              lite: <output folder>/colab_data_lite/ (under 1 GB), kept
The output folder is My Drive/CS3_workflow when Drive is connected, else the workflow folder on the Colab machine.
"""
import os
import resource
import shutil
import subprocess
import sys
from pathlib import Path

TOOLS = Path("/content/tools")
RMATS_ENV = Path("/content/rmats_env")
CORE = ["numpy", "scipy", "pandas", "matplotlib"]        # Colab has these loaded; pip must not replace them


def _sh(cmd, dry_run=False):
    """Prints the shell command `cmd` and runs it (only prints it if `dry_run`)."""
    print("$", cmd, flush=True)
    if not dry_run:
        subprocess.run(cmd, shell=True, check=True)


def _core_versions():
    """Returns the installed versions of numpy, scipy, pandas and matplotlib, so pip keeps Colab's versions."""
    import importlib.metadata as md
    return {p: md.version(p) for p in CORE}


def install_tools(dry_run=False):
    """samtools, bedtools and bedGraphToBigWig, plus pysam, pyBigWig, pybedtools, pyliftover and rbp-maps (~3 min)."""
    if (TOOLS / ".done").exists():
        return
    TOOLS.mkdir(parents=True, exist_ok=True)
    _sh("apt-get -qq update > /dev/null && apt-get -qq install -y samtools bedtools > /dev/null", dry_run)
    _sh(f"wget -q -O {TOOLS}/bedGraphToBigWig https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig"
        f" && chmod +x {TOOLS}/bedGraphToBigWig", dry_run)
    pins = TOOLS / "constraints.txt"
    pins.write_text("".join(f"{p}=={v}\n" for p, v in _core_versions().items()))
    _sh(f"pip install -q -c {pins} pysam pyBigWig pybedtools pyliftover seaborn tqdm", dry_run)
    # rbp-maps without its own version pins (it asks for newer numpy/scipy than Colab ships; Colab's work)
    _sh("pip install -q --no-deps --ignore-requires-python git+https://github.com/YeoLab/rbp-maps.git", dry_run)
    if not dry_run:
        (TOOLS / ".done").write_text("ok\n")


def install_rmats(version="4.3.0", dry_run=False):
    """rMATS-turbo from bioconda, in its own environment so Colab's Python is left alone (~5 min)."""
    if not (RMATS_ENV / "bin" / "rmats.py").exists():
        _sh(f"cd {TOOLS} && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba", dry_run)
        _sh(f"{TOOLS}/bin/micromamba create -y -q -r /content/micromamba -p {RMATS_ENV} -c conda-forge -c bioconda rmats={version}", dry_run)
    script = os.path.realpath(RMATS_ENV / "bin" / "rmats.py")
    os.environ["RMATS_CMD"] = f"{RMATS_ENV}/bin/python {script}"       # notebook 05 runs rMATS with this


def setup(workflow_folder, lite=False, data_on_drive=None, rmats=False, output_folder=None, dry_run=False):
    """Prepare this Colab session: install tools, set the two folders, raise the open-file limit.

    workflow_folder: the folder that holds these notebooks, workflow_utils.py and figure.py (the clone of the
        repository, or a copy on Google Drive).
    lite: use the small data in lite_data/ (fits free Colab) instead of downloading everything.
    data_on_drive: keep downloads and computed files on Drive instead of the Colab machine, so they survive the
        session. Default: yes in lite mode (under 1 GB), no otherwise (the full data needs about 50 GB).
    rmats: also install rMATS (notebook 05).
    output_folder: where results, the figure and (lite) the data go; default: workflow_folder.
    """
    workflow_folder = Path(workflow_folder).resolve()
    output_folder = Path(output_folder or workflow_folder).resolve()
    output_folder.mkdir(parents=True, exist_ok=True)
    assert (workflow_folder / "workflow_utils.py").exists(), f"workflow_utils.py not found in {workflow_folder}"
    if lite:
        assert (workflow_folder / "lite_data" / "MANIFEST.tsv").exists(), (
            f"lite_data/ not found in {workflow_folder}: copy it there, or set LITE = False")
        os.environ["CS3_LITE"] = "1"
        os.environ["CS3_LITE_DIR"] = str(workflow_folder / "lite_data")
    else:
        os.environ.pop("CS3_LITE", None)
    if data_on_drive is None:
        data_on_drive = lite
    os.chdir(workflow_folder)
    if str(workflow_folder) not in sys.path:
        sys.path.insert(0, str(workflow_folder))
    install_tools(dry_run)
    if rmats:
        install_rmats(dry_run=dry_run)
    os.environ["PATH"] = f"{TOOLS}:{os.environ['PATH']}"
    os.environ["CASE_STUDY_ROOT"] = str(output_folder / ("colab_output_lite" if lite else "colab_output"))
    data = output_folder / ("colab_data_lite" if lite else "colab_data") if data_on_drive else Path("/content/cs3_data")
    data.mkdir(parents=True, exist_ok=True)
    os.environ["CASE_STUDY_DATA"] = str(data)
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)          # rbp-maps keeps many files open
    resource.setrlimit(resource.RLIMIT_NOFILE, (min(10240, hard), hard))
    free = shutil.disk_usage(data).free / 1e9
    print("lite mode (small data from lite_data/)" if lite else "full data")
    print(f"results and figure: {os.environ['CASE_STUDY_ROOT']}")
    print(f"downloads: {data} ({free:.0f} GB free)")
