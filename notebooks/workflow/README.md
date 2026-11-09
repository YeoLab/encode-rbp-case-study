# Case Study 3 workflow notebooks

These seven notebooks follow figure panel a, one notebook per box. Together they rebuild the Case Study 3 figure: QKI regulates two PICALM cassette exons, the 150-nt exon (ENSE00000989276.1) and the 24-nt exon (ENSE00001416869.1), in opposite directions.

| Notebook | Box in panel a | What it does | Main output |
| --- | --- | --- | --- |
| 00_colab_run_all | | On Google Colab, runs 01 to 07 in one session and saves each finished notebook to Drive | `colab_output/executed/` |
| 01_gene_search | Gene search (RBP-ARK) | Screenshots of rbp-ark.com's gene search for PICALM, then the same search through the API | `results/workflow/01_gene_search_PICALM.tsv` |
| 02_experiment_detail | Experiment detail (RBP-ARK) | QC, track hubs and windows files of the three QKI eCLIP experiments | `results/workflow/02_qki_windows_PICALM.tsv` |
| 03_download | Download (RBP-ARK) | BDBag, API file list and Globus; downloads the eCLIP BAMs the splice maps need (about 4 GB) | `results/workflow/03_eclip_bams.tsv` |
| 04_shrna_rnaseq | shRNA RNA-seq (ENCODE) | The eight QKI knockdown and control RNA-seq BAMs (about 40 GB) and the QKI knockdown level | `results/workflow/04_knockdown_level.tsv` |
| 05_differential_splicing | Differential splicing (ENCODE) | rMATS, then the three exon groups lifted from hg19 to GRCh38 | event lists in `<data>/rbp_maps/event_lists_rerun/` |
| 06_dpsi_per_exon | ΔPSI per exon (ENCODE) | Junction reads, PSI, ΔPSI and FDR of the two PICALM exons | `results/workflow/06_picalm_knockdown.json` |
| 07_integrate | Integrate | Isoforms, rbp-maps splice maps and the figure | `figures/workflow/CS3_PICALM_QKI_workflow.png` and `.pdf` |

Run them in order, 01 to 07. On Google Colab, `00_colab_run_all.ipynb` runs all seven for you.

## Run on Google Colab

You need a Google account. Nothing needs to be installed on your computer, and the default CPU runtime is enough (no GPU).

### 1. Open the notebook

Open [00_colab_run_all.ipynb in Colab](https://colab.research.google.com/github/yeolab/encode-rbp-case-study/blob/main/notebooks/workflow/00_colab_run_all.ipynb). Nothing needs to be uploaded: the first code cell copies this repository, with the small data in `lite_data/`, onto the Colab machine.

### 2. Why this notebook

Each Colab notebook runs on its own machine, and that machine's files are deleted when the session ends. Notebooks 05 to 07 need the large files that 03 and 04 download, so they have to run on the same machine. `00_colab_run_all` runs notebooks 01 to 07 one after another on its machine.

### 3. Choose lite or full

The first code cell of every notebook has `LITE = True`. That is the free Colab version: it uses the small files in `lite_data/` instead of downloading about 50 GB, and the whole run takes about an hour. Set `LITE = False` to download and process everything (see "Lite and full" below).

### 4. Run it

Choose Runtime > Run all.

- The first code cell copies the repository and asks for permission to connect your Google Drive. Allow it, so the results and the figure are saved in `My Drive/CS3_workflow/`. To keep them only on the Colab machine instead, set `SAVE_TO_DRIVE = False`.
- It then installs samtools, bedtools, bedGraphToBigWig, rbp-maps, rMATS and a few Python packages (about 10 minutes) and prints how much disk space the Colab machine has.
- The next cell runs the seven notebooks and prints a line as each one starts and finishes.

The lite run takes about an hour; the full run takes 5 to 8 hours. Keep the browser tab open and your computer awake: Colab stops a session that looks idle, and free sessions last at most 12 hours. A paid Colab plan gives longer sessions that keep running when you close the tab.

`REP1_ONLY` computes only replicate 1 of the splice maps, the one the figure shows. It is on in lite mode. In full mode it saves about 40 minutes.

### 5. Look at the results

In lite mode the folder is `colab_output_lite/` instead of `colab_output/`.

| What | Where in `My Drive/CS3_workflow/` |
| --- | --- |
| The figure | `colab_output/figures/workflow/CS3_PICALM_QKI_workflow.png` and `.pdf` |
| Each notebook, run, with its tables and plots | `colab_output/executed/01_gene_search.ipynb` to `07_integrate.ipynb` (open them from Drive with Colaboratory to read them) |
| Results tables | `colab_output/results/workflow/` |

In full mode, downloads, rMATS runs and splice maps stay on the Colab machine (`/content/cs3_data/`) and are deleted when the session ends. In lite mode they are small, so they go to `colab_data_lite/` on Drive and are kept.

### Lite and full

The figure needs only a small part of the 50 GB of BAM files. The lite version in `lite_data/` (about 90 MB) keeps just that part:

| Data | Full | Lite (`lite_data/`) |
| --- | --- | --- |
| eCLIP BAMs (notebooks 03, 07) | Seven BAMs, about 4 GB | About 1% of the reads: read 2 within 400 nt of every splice site of the exons in the splice maps, plus the PICALM gene. rbp-maps reads 350 nt around each site, so nothing it uses is lost. `library_sizes.tsv` keeps each full BAM's read count, so the signal is scaled as if it came from the full BAM. |
| shRNA RNA-seq BAMs (04, 05) | Eight BAMs, about 40 GB | Only the reads on PICALM (hg19, ±2 kb) |
| rMATS (05) | Genome-wide, 2 to 4 h | PICALM only, with PICALM's gene models, about a minute |
| Exon groups for the splice maps (05, 07) | From the genome-wide rMATS run | From the full run's tables, shipped in `lite_data/rmats_full/` |

What stays the same: PICALM's junction reads, PSI and ΔPSI (notebook 05 prints the lite and full rows side by side), the exon groups, and the splice maps. What differs: rMATS's P values and FDR for PICALM, because the lite run tests a handful of exons instead of tens of thousands. Notebook 06 and the figure use the full run's FDR from `lite_data/rmats_full/`.

Notebooks 01 and 02 query rbp-ark.com directly and are the same in both modes.

#### Making lite_data/ (maintainers)

`lite_data/` is built from the full data, once, on a computer that has run notebooks 03, 04 and 05 in full:

```
conda activate encode-rbp-case-study
python make_lite_data.py ../../data      # or the folder CASE_STUDY_DATA points to
```

It writes `lite_data/` next to the script, with `MANIFEST.tsv` (sizes and MD5 checksums). Commit it with the rest of the repository.

### If the run stops partway

The cell prints the error and which notebook it came from. In the same session, fix the problem, set `START` to that notebook's number (for example `START = "05"`) and run the cell again; finished downloads and computations on the machine are reused. In lite mode everything is on Drive, so this also works in a new session. In full mode, if the session has ended, the machine's files are gone: set `START = "03"` (01 and 02 only save small tables to Drive) and run again.

Common problems:

- "workflow_utils.py not found" or "lite_data/ not found": the copy of the repository on the Colab machine is incomplete. Delete it (`!rm -rf /content/encode-rbp-case-study`) and run the first cell again.
- Out of disk space (full mode): notebooks 03 to 07 need about 50 GB on the Colab machine, and the first cell prints how much is free. Use a runtime with more disk, or keep the data on Drive (below).

### Running notebooks one at a time

You can also open and run each notebook on its own (in Colab: File > Open notebook > GitHub, `yeolab/encode-rbp-case-study`): every notebook starts with the same Colab cell. Notebooks 01 and 02 work this way without changes, and so do all of them in lite mode, which keeps its files on Drive. In full mode, from 03 on, the notebooks pass large files to each other, so they must be kept on Drive rather than on each notebook's machine: in each notebook's first code cell, add `data_on_drive=True` to the `colab_setup.setup(...)` call, for example `colab_setup.setup(WORKFLOW, lite=LITE, output_folder=OUTPUT, data_on_drive=True)`. This needs about 50 GB of free Drive space. Everything finished is then kept between sessions, so you can stop and continue another day.

### Time and disk, notebook by notebook (full mode)

Rough times on Colab's standard CPU runtime (two processor cores):

| Notebook | Downloads | Rough time | Needs from earlier notebooks |
| --- | --- | --- | --- |
| 01 | under 1 MB | 1 min | nothing |
| 02 | about 30 MB | 2 min | nothing |
| 03 | about 4 GB of eCLIP BAMs | 15 to 30 min | nothing |
| 04 | about 40 GB of RNA-seq BAMs | 1 to 2 h | nothing |
| 05 | 1.8 GB of gene annotation | 2 to 4 h (rMATS, twice) | 04's BAMs |
| 06 | none | 1 min | 05's rMATS results |
| 07 | about 50 MB | 1 to 1.5 h (four splice maps) | 03's BAMs, 05's event lists |

## Run on your own computer

Use the `encode-rbp-case-study` environment from the case study's `notebooks/environment.yml` (rMATS 4.3.0, rbp-maps, samtools, bedtools, bedGraphToBigWig, pysam, pyliftover). Open the notebooks from this folder so that `workflow_utils.py` and `figure.py` are found. The Colab cell at the top of each notebook does nothing outside Colab.

Two folders:

- The case study folder (`CASE_STUDY_ROOT`, by default two levels up from here) gets the small results and the figure.
- The data folder (`CASE_STUDY_DATA`, by default `data/` in the case study folder) holds downloads, rMATS runs and splice maps. These notebooks use the same places as notebooks 00 to 05, so files already downloaded or computed there are reused.

To run all seven from a terminal: `bash run_workflow.sh`.

## Rerunning

Every download is cached, and large files are checked against their checksums once. A finished rMATS run is reused when its annotation, `--libType` and rMATS version match. A finished splice map is reused when all its outputs are newer than its event lists; event lists are only rewritten when their content changes.

## Files

- `workflow_utils.py`: the functions the notebooks call, each with a docstring.
- `figure.py`: `draw_figure()`, which draws panels a to e.
- `colab_setup.py`: installs the tools and sets the folders on Google Colab.
- `lite_data/` and `make_lite_data.py`: the small data for free Colab, and the script that makes it.
- `build_notebooks.py` and `screenshots/`: rebuild the notebooks (`python build_notebooks.py .`). The screenshots are already embedded in the notebooks.
