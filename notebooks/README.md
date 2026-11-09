# Case Study 3 notebooks: QKI regulates two PICALM exons in opposite directions

These notebooks rebuild the Case Study 3 figure from public data (download the knockdown BAMs → rerun rMATS → sashimi plots → eCLIP BAMs from RBP-ARK → rbp-maps). Run them in order, from this `notebooks/` folder. The per-box notebooks in `workflow/` reproduce the same figure, one notebook per box of panel a.

| # | Notebook | What it does | Main outputs | Time / disk |
| --- | --- | --- | --- | --- |
| 00 | `00_screen_and_locus.ipynb` | Panels a, b: screens 18 RBPs for ENCORE/ENCORE2 binding next to exons that respond to knockdown in both cell lines; QKI windows, RNA and transcripts at PICALM's 150-nt and 24-nt exons | `results/screen_candidates.tsv`, `figures/data/picalm_locus.json` | ~0.4 GB; 10-20 min |
| 01 | `01_download_shrna_bams.ipynb` | Reads from the ENCODE API which BAMs ENCODE's QKI rMATS results were built from, then downloads and indexes those 8 BAMs (hg19). Also fetches ENCODE's rMATS releases and both DESeq2 tables per knockdown | `<data>/downloads/shrna_bams_hg19/` (+ `manifest.tsv`), `<data>/downloads/encode_reference/` | 38.7 GB |
| 02 | `02_run_rmats.ipynb` | Runs rMATS-turbo 4.3.0 (knockdown vs non-target control, GENCODE v40, `--libType fr-firststrand`) for QKI in K562 and HepG2, and compares with ENCODE's rMATS 3.2.1 release | `<data>/rmats/rerun/<knockdown>/SE.MATS.JC.txt`, `<data>/rmats/encode_release/<knockdown>/` | ~30-60 min per knockdown |
| 03 | `03_panel_c_sashimi.ipynb` | Panel c: sashimi plots of PICALM exon 13 and the 24-nt exon (both QKI knockdowns), PSI per replicate, ΔPSI, FDR, knockdown level | `figures/panel_c/`, `results/panel_c_stats_<source>.tsv` | minutes |
| 04 | `04_eclip_rbp_maps.ipynb` | Panels d, e: downloads the 7 eCLIP BAMs of the two QKI datasets in the figure's splice maps from RBP-ARK, builds the event lists (one event per exon with rbp-maps' `subset_rmats_junctioncountonly.py`, then lifted to GRCh38), makes the splice maps with rbp-maps (both datasets against the HepG2 QKI knockdown), and extracts PICALM's own traces | `<data>/rbp_maps/`, `figures/panel_de/`, `results/rbp_maps_*` | 3.7 GB; ~10 min per map (4 maps) with the permutation test |
| 05 | `05_window_maps_optional.ipynb` | Optional: coarser splice maps from RBP-ARK's windows table, no BAMs | `results/window_map_*` | minutes |

Notebooks 03 and 04 can use either the rMATS rerun from notebook 02 (default) or ENCODE's released rMATS output, which the draft figure used:

```bash
SE_SOURCE=encode jupyter lab
```

## Annotation

All analyses use **GENCODE v40**. The knockdown BAMs are hg19, so rMATS (notebook 02) uses GENCODE's GRCh37 mapping of v40 (`gencode.v40lift37.annotation.gtf`); transcripts in panel b come from the GRCh38 v40 GTF (notebook 00). ENCODE's released rMATS tables (`SE_SOURCE=encode`, and the screen in notebook 00) were computed by ENCODE with GENCODE v19. RBP-ARK's window annotations (feature types) come from Skipper.

## Folders

Small results and figures go in this case study folder. Downloads and other large files (about 45 GB) go in a data folder, by default `data/` in the case study folder (set `CASE_STUDY_DATA` to keep them elsewhere, such as an external drive):

| Path in the data folder | Contents |
| --- | --- |
| `downloads/screen/` | RBP-ARK catalog and windows table, ENCODE rMATS tars for the screen, liftover chain, GENCODE v40 GTF (notebooks 00, 05) |
| `downloads/shrna_bams_hg19/`, `downloads/encode_reference/` | knockdown BAMs, rMATS releases, DESeq2 tables, GENCODE v40 GRCh37 mapping (01-03) |
| `rmats/` | rMATS reruns and extracted releases (02) |
| `rbp_ark_eclip_skipper2/full/` | the eCLIP BAMs (04; `../downloads/get_data.sh` fills the same folder) |
| `rbp_maps/` | event lists, bigWigs, maps (04) |

## Setup

```bash
conda env create -f environment.yml
conda activate encode-rbp-case-study
python -m ipykernel install --user --name encode-rbp-case-study
jupyter lab
```

The environment provides rMATS-turbo, rmats2sashimiplot, rbp-maps, samtools, bedtools, bedGraphToBigWig, pysam, pyBigWig and pyliftover. Network access is needed for encodeproject.org, the RBP-ARK data API and its Globus data endpoint, api.genome.ucsc.edu, ftp.ebi.ac.uk (GENCODE v40) and hgdownload.soe.ucsc.edu (liftover chain).

Environment variables (all optional):

| Variable | Default | Used by |
| --- | --- | --- |
| `CASE_STUDY_ROOT` | `..` (the case study folder) | all |
| `CASE_STUDY_DATA` | `../data` (`data/` in the case study folder) | all |
| `RBPARK_API` | `https://rbp-ark.com` (the RBP-ARK production site and public data API) | 00, 04, 05 |
| `SE_SOURCE` | `rerun` | 03, 04 |
| `RMATS_THREADS` | all cores | 02 |
| `RBPMAPS_SIGTEST` | `permutation` (`zscore` is much faster for a first look) | 04 |

## How the numbers should compare with the draft figure

- **Panel c.** With `SE_SOURCE=encode`, exon 13 ΔPSI is +0.22 (K562) and +0.48 (HepG2); 24-nt exon −0.12 and −0.53. The rMATS-turbo rerun should agree to within a few reads.
- **Knockdown level.** From ENCODE's DESeq2 tables without batch normalization: QKI log2FC −0.85 (K562) and −3.36 (HepG2). The batch-normalized tables give −1.21 and −6.06. Notebook 03 prints both.
- **Event lists.** Rules (the same as the earlier manual run): included / excluded = FDR < 0.1, P < 0.05, ΔPSI > 0.05 / < −0.05; background = P > 0.5 with mean PSI strictly between 0.05 and 0.95 in both knockdown and control; no read minimum. Each set then keeps one event per exon by the rule of rbp-maps' `subset_rmats_junctioncountonly.py` (`cs3_utils.one_event_per_exon()`; a tie goes to the event listed first, where rbp-maps' result depends on the computer), and is lifted to GRCh38 (each exon whole, same chromosome and strand, unchanged length).
  - ENCODE release (hg19 → one per exon, included / excluded / background): QKI K562 18→15 / 15→15 / 5,412→3,737; QKI HepG2 390→298 / 130→107 / 4,693→3,334.
  - Rerun: counts are printed by notebook 04 (`results/rbp_maps_event_counts_rerun.tsv`).
  - Checked against rbp-maps' script on the rerun tables: the same events in every set except where two forms of an exon tie (2 to 12 per set), which rbp-maps itself resolves differently on macOS and Linux.
- **Panels d, e.** RBP-ARK's windows put QKI 0-100 nt upstream of exon 13 and 0-100 nt downstream of the 24-nt exon; the PICALM exons' own traces should peak there (exon 13: 33 nt upstream with the K562 eCLIP, 11 nt with the HepG2 eCLIP).

## Tested

All six notebooks were executed end to end with nbclient on a synthetic dataset with the same file layout, served through stand-ins for `requests` and `curl` (ENCODE file objects with `derived_from` lineage and MD5s, rMATS release tars, DESeq2 tables, RBP-ARK API file lists with SHA-256s, UCSC JSON, a liftover chain), with rMATS-turbo 4.4.0, rmats2sashimiplot, rbp-maps and Python 3.12. The splice maps recovered the binding positions planted in the synthetic eCLIP reads. The event-list counts above were computed from ENCODE's real releases. The notebooks were then run on the full ENCODE and RBP-ARK data (Nov 2026), which produced the figure's inputs. Their saved outputs are not included; rerunning reuses any downloads and computed files already in the data folder.
