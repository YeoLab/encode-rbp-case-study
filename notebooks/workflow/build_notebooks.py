"""Build the seven Case Study 3 workflow notebooks (one per box of figure panel a).

    python build_notebooks.py [out_dir]

Screenshots in screenshots/ are embedded in the notebooks as attachments.
"""
import base64
import sys
from pathlib import Path

import nbformat

HERE = Path(__file__).parent
SHOTS = HERE / "screenshots"

SETUP = '''import json
import sys
from pathlib import Path

import pandas as pd
from IPython.display import Image, Markdown, display

sys.path.insert(0, ".")
import workflow_utils as U

ROOT, DATA = U.folders()
RESULTS = ROOT / "results" / "workflow"
API_CACHE = DATA / "downloads" / "rbpark_api"
print("case study folder:", ROOT)
print("data folder:", DATA)'''


COLAB_MD = """On Google Colab, run the next cell first in every session. It copies this repository onto the Colab machine and, with `SAVE_TO_DRIVE = True`, saves results and the figure to `My Drive/CS3_workflow` (Colab asks for permission to connect your Drive). `LITE = True` uses a small subset of the data that fits free Colab. On your own computer the cell does nothing."""

COLAB = '''import sys

if "google.colab" in sys.modules:
    import os
    import subprocess
    REPO = "https://github.com/yeolab/encode-rbp-case-study"
    CLONE = "/content/encode-rbp-case-study"
    LITE = True    # True: small subset of the data, fits free Colab (about 1 h). False: all the data (5-8 h, ~50 GB)
    SAVE_TO_DRIVE = True    # results and the figure in My Drive/CS3_workflow; False: on this Colab machine only
    if not os.path.isdir(CLONE):
        subprocess.run(["git", "clone", "-q", "--depth", "1", REPO, CLONE], check=True)
    WORKFLOW = CLONE + "/notebooks/workflow"
    OUTPUT = None
    if SAVE_TO_DRIVE:
        from google.colab import drive
        drive.mount("/content/drive")
        OUTPUT = "/content/drive/MyDrive/CS3_workflow"
    sys.path.insert(0, WORKFLOW)
    import colab_setup
    colab_setup.setup(WORKFLOW, lite=LITE, output_folder=OUTPUT{extra})'''


def shot(name, alt):
    """Returns a screenshot cell: the image screenshots/`name` with alt text `alt`."""
    return ("shot", name, alt)


NB = {}

# ------------------------------------------------------------------ 01 gene search
NB["01_gene_search"] = [
    ("md", """# 01 · Gene search (RBP-ARK)

Which RNA-binding proteins (RBPs) bind PICALM, and where? RBP-ARK's gene search lists their reproducible enriched windows. A window is a stretch of about 100 nt where the eCLIP signal of an RBP is higher than in its size-matched input, in both replicates.

This notebook shows the search in the browser, then runs the same search in code. Notebook 02 then looks at the QKI experiments one by one."""),
    ("md", """## In the browser

1. Open [rbp-ark.com/genesearch](https://rbp-ark.com/genesearch/), type PICALM and keep the minimum enrichment at ≥ 3 (mean log2 fold over input). The By RBP tab ranks the 120 experiments with a window on PICALM by their number of windows."""),
    shot("01a_gene_search_by_rbp.jpg", "Gene search results for PICALM, By RBP tab"),
    ("md", """2. Open the All windows tab and type QKI in the filter box under Experiment. Both K562 experiments have a window at 85,990,399-85,990,499, next to the 150-nt exon, and all three have one at 85,977,969-85,978,068, next to the 24-nt exon."""),
    shot("01b_gene_search_qki_windows.jpg", "All windows tab filtered to QKI"),
    ("md", """3. Scroll down to the coordinate view, pick "Selected gene, all isoforms", and drag across the plot to zoom in on 85,976,000-85,992,500. The top rows are PICALM's transcripts and the bottom rows the QKI windows. PICALM is on the minus strand, so it is read from right to left here."""),
    shot("01c_coordinate_view_qki.jpg", "Coordinate view zoomed on the two PICALM exons"),
    ("md", """## The same search in code

The gene search page reads RBP-ARK's public API. `gene_search()` asks `/api/genes/PICALM/experiments/?min_enrichment=3` and returns one row per window and experiment."""),
    ("code", SETUP),
    ("code", """hits = U.gene_search(U.GENE, API_CACHE, min_enrichment=3)
print(hits["subject"].nunique(), "experiments,", len(hits), "windows on", U.GENE)
hits.groupby(["subject", "rbp"]).size().sort_values(ascending=False).head(10).rename("windows")"""),
    ("md", """The QKI windows, with each window's position next to the two exons. Positions count like intron nucleotides: −1 is the last nucleotide before an exon's 3′ splice site, +1 the first one after its 5′ splice site."""),
    ("code", """qki = hits[hits["rbp"] == U.RBP].copy()
for exon, info in U.EXONS.items():
    qki["nt from " + info["name"]] = [U.relative_position(s, e, exon) for s, e in zip(qki["start"], qki["end"])]
qki.to_csv(RESULTS / "01_gene_search_PICALM.tsv", sep="\\t", index=False)
qki[["subject", "start", "end", "feature_type", "log2_enrichment", "q_value", "nt from 150-nt exon", "nt from 24-nt exon"]]"""),
    ("md", """Two windows sit right next to the exons: one ending 1 nt before the 150-nt exon (both K562 experiments) and one starting 2 nt after the 24-nt exon (all three). The HepG2 experiment also has a window before the 150-nt exon, but below the search's cutoff of 3 (notebook 02). Notebook 07 checks whether these positions match what QKI does across the genome."""),
]

# ------------------------------------------------------------------ 02 experiment detail
NB["02_experiment_detail"] = [
    ("md", """# 02 · Experiment detail (RBP-ARK)

Every eCLIP experiment has its own page on RBP-ARK, with quality control (QC), a UCSC track hub and the file of windows. QKI has three experiments that pass QC: one from ENCORE2 in K562 cells and two from ENCODE3, in K562 and in HepG2 cells.

This notebook collects the same information in code and saves each experiment's windows near PICALM, which notebook 07 draws in the figure."""),
    ("md", """## In the browser

1. On the [Experiments](https://rbp-ark.com/search/?q=QKI) page, search for QKI. All three experiments pass QC."""),
    shot("02a_experiments_qki.jpg", "Experiments page searched for QKI"),
    ("md", """2. Open one, here the ENCORE2 experiment (accession 5007). The read table shows how many reads are left after trimming, mapping and removing PCR duplicates. The buttons at the top download a BDBag (notebook 03) or start a Globus transfer."""),
    shot("02b_experiment_detail.jpg", "Experiment page of QKI_K562_5007"),
    ("md", """3. Further down are the track hub address and the windows file. To see the signal in the UCSC Genome Browser, paste the address into My Data > Track Hubs (or use the links the code below prints)."""),
    shot("02c_trackhub_and_windows.jpg", "Track hub address and windows section"),
    ("md", "## In code"),
    ("code", SETUP),
    ("code", """experiments = U.experiments_table(API_CACHE)
experiments[experiments["rbp"] == U.RBP]"""),
    ("md", """Reads per sample, from the same QC summary the page shows. IP is the immunoprecipitated sample; IN is its size-matched input."""),
    ("code", """stats = pd.concat({subject: U.read_stats(subject, API_CACHE) for subject in U.ECLIP}, names=["experiment"])
stats"""),
    ("code", """region = "{}:{}-{}".format(U.LOCUS[0], U.LOCUS[2], U.LOCUS[3])
for subject in U.ECLIP:
    hub = U.experiment_record(subject, API_CACHE, "trackhub_url")["trackhub_url"]
    link = f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubUrl={hub}&position={region}"
    display(Markdown(f"{subject}: [open the track hub at PICALM in UCSC]({link})"))"""),
    ("md", """## Windows near PICALM

`reproducible_windows()` downloads an experiment's windows file (a few MB) once. `log2_enrichment` is the mean log2 fold over input and `q_value` the adjusted p-value; these are the numbers the gene search shows. The figure shades windows by enrichment and outlines those at log2 ≥ 3."""),
    ("code", """parts = []
for subject in U.ECLIP:
    near = U.windows_in(U.reproducible_windows(subject, API_CACHE), *U.LOCUS)
    parts.append(near.assign(subject=subject)[["subject", "chr", "start", "end", "strand", "feature_type_top",
                                                "log2_enrichment", "q_value"]])
windows = pd.concat(parts, ignore_index=True)
for exon, info in U.EXONS.items():
    windows["nt from " + info["name"]] = [U.relative_position(s, e, exon) for s, e in zip(windows["start"], windows["end"])]
windows.to_csv(RESULTS / "02_qki_windows_PICALM.tsv", sep="\\t", index=False)
print(len(windows), "QKI windows between", U.LOCUS[2], "and", U.LOCUS[3])
windows[windows[["nt from 150-nt exon", "nt from 24-nt exon"]].abs().min(axis=1) <= 100]"""),
    ("md", """Expected: the window before the 150-nt exon has log2 enrichment 4.08 (ENCORE2), 3.31 (ENCODE3 K562) and 2.40 (HepG2); the window after the 24-nt exon has 4.59, 3.40 and 4.28."""),
]

# ------------------------------------------------------------------ 03 download
NB["03_download"] = [
    ("md", """# 03 · Download (RBP-ARK)

The splice maps in notebook 07 need the eCLIP reads (BAM files), not only the windows. RBP-ARK offers three ways to get files:

| Way | Good for | How |
| --- | --- | --- |
| BDBag | one experiment | the Download BDBag button on its page: a small archive that lists every file with its address, size and checksum |
| API file list | a few files, from a script | what this notebook does |
| Globus | many experiments or whole folders | log in with Globus and transfer to your computer or cluster |

Please download only what you need. The public server is shared, and bulk downloads have taken it offline before."""),
    shot("03a_download_options.jpg", "Download BDBag and Globus transfer on an experiment page"),
    ("code", SETUP + '''
import subprocess
import tarfile'''),
    ("md", """## A BDBag

A BDBag is a few kB. Its `fetch.txt` lists each file's address, size and path."""),
    ("code", """bag = U.fetch(f"{U.RBPARK}/search/detail/QKI_K562_5007/bdbag/", API_CACHE / "QKI_K562_5007.bdbag.tar.gz")
with tarfile.open(bag) as tar:
    fetch_txt = tar.extractfile(next(m for m in tar.getmembers() if m.name.endswith("fetch.txt"))).read().decode()
listed = pd.DataFrame([line.split("\\t") for line in fetch_txt.splitlines()], columns=["url", "size", "path"])
listed["size"] = listed["size"].astype("int64")
print(len(listed), "files, {:.1f} GB in all".format(listed["size"].sum() / 1e9))
listed[listed["path"].str.contains("/bams/dedup/")]"""),
    ("md", """To download a whole bag from a terminal (about 6 GB for this experiment):

```bash
pip install bdbag
tar xzf QKI_K562_5007.bdbag.tar.gz
bdbag --resolve-fetch all QKI_K562_5007
```

## Globus, for many files

Sign in at [app.globus.org](https://app.globus.org), open the RBP-ARK collection (`9cd9be17-cbed-4b34-ab76-704b2fd6ac7b`) and transfer folders to your own endpoint, for example your laptop with [Globus Connect Personal](https://www.globus.org/globus-connect-personal). The same from a terminal:

```bash
pip install globus-cli
globus login
RBPARK=9cd9be17-cbed-4b34-ab76-704b2fd6ac7b
OUT=/encore2/eclip/Skipper/2.0.0/v2/output
globus ls "$RBPARK:$OUT/secondary_results/bams/dedup/genome/" | grep QKI
globus transfer "$RBPARK:$OUT/secondary_results/bams/dedup/genome/QKI_K562_5007_IP_1.genome.Aligned.sort.dedup.bam" \\
    "$(globus endpoint local-id):$HOME/rbpark/QKI_K562_5007_IP_1.genome.Aligned.sort.dedup.bam"
```

## The BAMs for the splice maps

Each splice map compares an IP sample with its size-matched input. `eclip_bams()` finds the deduplicated BAMs of both replicates in each experiment's file list. `download_large()` resumes a partial download and checks the SHA-256 checksum, so a file already on disk is not downloaded again.

In lite mode (free Colab) nothing is downloaded here. The lite BAMs in `lite_data/eclip/` hold only the reads the splice maps use: read 2 within 400 nt of every splice site of the exons in the maps, a small part of each file. Next to them, `library_sizes.tsv` keeps each full BAM's read count, so notebook 07 scales the signal as if it came from the full BAM."""),
    ("code", """bams = U.eclip_bams(API_CACHE)
print("{} files, {:.1f} GB".format(len(bams), bams["size"].sum() / 1e9))
bams[["subject", "sample", "size", "collection_path"]]"""),
    ("code", """if U.LITE:
    lite = U.lite_eclip_bams()
    bams = bams.merge(lite[["subject", "sample", "bam", "mapped_F260", "mapped_records", "kept_reads"]], on=["subject", "sample"])
    print("lite mode: using", len(bams), "BAMs from", U.LITE_DIR / "eclip")
else:
    ECLIP_DIR = DATA / "rbp_ark_eclip_skipper2" / "full"
    local = []
    for _, row in bams.iterrows():
        out = ECLIP_DIR / Path(row["collection_path"]).name
        print("checking" if out.exists() else "downloading", out.name, flush=True)
        U.download_large(row["url"], out, row["size"], row["sha256"])
        if not Path(str(out) + ".bai").exists():
            subprocess.check_call(["samtools", "index", str(out)])
        local.append(str(out))
    bams["bam"] = local
bams.to_csv(RESULTS / "03_eclip_bams.tsv", sep="\\t", index=False)
bams.drop(columns=["url", "collection_path", "sha256"])"""),
]

# ------------------------------------------------------------------ 04 shRNA RNA-seq
NB["04_shrna_rnaseq"] = [
    ("md", """# 04 · shRNA RNA-seq (ENCODE)

To see what QKI does to splicing, ENCODE knocked QKI down with a short hairpin RNA (shRNA) in K562 and HepG2 cells and sequenced their RNA. Control cells got a shRNA that targets nothing. Each condition has two replicates, so there are eight RNA-seq files (BAMs, aligned to the hg19 genome)."""),
    shot("04a_encode_qki_hepg2_knockdown.jpg", "ENCODE page of the QKI knockdown in HepG2"),
    ("md", """## Which BAMs

ENCODE's own splicing analysis (rMATS) records which BAMs it used. `knockdown_bams()` reads that list from the ENCODE API, so this notebook starts from the same files."""),
    ("code", SETUP + '''
import subprocess

ENCODE_CACHE = DATA / "downloads" / "encode_api"
BAM_DIR = DATA / "downloads" / "shrna_bams_hg19"
REF_DIR = DATA / "downloads" / "encode_reference"'''),
    ("code", """lineage = pd.concat([U.knockdown_bams(key, ENCODE_CACHE) for key in U.KNOCKDOWNS], ignore_index=True)
print("{:.0f} GB in all".format(lineage["file_size"].sum() / 1e9))
lineage[["knockdown", "condition", "rep", "experiment", "accession", "assembly", "read_length", "file_size"]]"""),
    ("md", """Expected: K562 knockdown ENCFF280ZIB and ENCFF708KBX, control ENCFF092XAP and ENCFF265ZZB; HepG2 knockdown ENCFF610CID and ENCFF262GAC, control ENCFF755GKJ and ENCFF571ARP.

The next cell downloads them (about 40 GB) and indexes them. It checks each file's MD5 checksum, which takes about a minute per GB the first time. In lite mode it uses `lite_data/shrna/` instead: the same eight BAMs cut down to the reads on PICALM, a few MB each."""),
    ("code", """paths = []
for _, row in lineage.iterrows():
    if U.LITE:
        paths.append(str(U.lite_shrna_bam(row["accession"])))
        continue
    out = BAM_DIR / (row["accession"] + ".bam")
    print("checking" if out.exists() else "downloading", row["knockdown"], row["condition"], "rep", row["rep"], out.name, flush=True)
    U.download_large(row["url"], out, row["file_size"], row["md5sum"], "md5")
    if not Path(str(out) + ".bai").exists():
        subprocess.check_call(["samtools", "index", str(out)])
    paths.append(str(out))
lineage["bam"] = paths
BAM_DIR.mkdir(parents=True, exist_ok=True)
lineage.to_csv(BAM_DIR / ("manifest_qki_lite.tsv" if U.LITE else "manifest_qki.tsv"), sep="\\t", index=False)
lineage[["knockdown", "condition", "rep", "bam"]]"""),
    ("md", """## How much QKI is left

ENCODE also measured gene expression in the same cells (DESeq2). The log2 fold change of QKI mRNA, knockdown vs control, shows how strong each knockdown is. ENCODE published two tables per knockdown, with and without batch normalization; this case study uses the one without."""),
    ("code", """rows = []
for key, kd in U.KNOCKDOWNS.items():
    for role in ("deseq2", "deseq2_batch", "rmats"):
        U.encode_file(kd[role], REF_DIR, ENCODE_CACHE)
    rows.append({"knockdown": key,
                 "QKI log2FC": U.deseq2_log2fc(REF_DIR / (kd["deseq2"] + ".tsv"), kd["gene_id"]),
                 "QKI log2FC, batch-normalized": U.deseq2_log2fc(REF_DIR / (kd["deseq2_batch"] + ".tsv"), kd["gene_id"])})
level = pd.DataFrame(rows)
level["QKI mRNA left"] = (2 ** level["QKI log2FC"]).map("{:.0%}".format)
level.to_csv(RESULTS / "04_knockdown_level.tsv", sep="\\t", index=False)
level"""),
    ("md", """QKI mRNA drops to about 55% of control in K562 (log2 −0.85) and to about 10% in HepG2 (log2 −3.36). The HepG2 knockdown is much stronger, so the later notebooks lean on HepG2."""),
]

# ------------------------------------------------------------------ 05 differential splicing
NB["05_differential_splicing"] = [
    ("md", """# 05 · Differential splicing (rMATS)

A cassette exon is an exon that cells sometimes include in the mRNA and sometimes skip. PSI (percent spliced in) is the fraction of transcripts that include it, from 0 to 1. rMATS counts the RNA-seq reads that cross exon junctions, both those that include the exon and those that skip it. From these it computes PSI per replicate and tests whether PSI differs between knockdown and control. ΔPSI is PSI after knockdown minus PSI in control.

This notebook runs rMATS on notebook 04's BAMs. It then sorts the exons into the three groups that the splice maps compare and moves them from hg19 to GRCh38, the genome version RBP-ARK uses."""),
    ("code", SETUP + '''
import os
import shutil
import subprocess

BAM_DIR = DATA / "downloads" / "shrna_bams_hg19"
REF_DIR = DATA / "downloads" / "encode_reference"
GENOME_CACHE = DATA / "downloads" / "genome"
DIRS = U.data_dirs(DATA)
EVENT_DIR = DIRS["events"]
THREADS = os.cpu_count() or 4'''),
    ("md", """## Run rMATS

rMATS needs gene models in the BAMs' genome version, so it uses GENCODE v40 mapped to hg19 (`gencode.v40lift37`). The reads are paired-end, 100 nt, from strand-specific libraries (`--libType fr-firststrand`). One run takes 30 to 60 minutes. A finished run with the same settings is reused.

In lite mode rMATS runs on the PICALM-only BAMs with PICALM's gene models, which takes a minute. It tests only PICALM's exons, so the splice maps further down use the genome-wide results of the full run, which come with `lite_data/`."""),
    ("code", """GTF = U.LITE_DIR / "gencode.v40lift37.PICALM.gtf" if U.LITE else REF_DIR / "gencode.v40lift37.annotation.gtf"
if not GTF.exists():
    U.fetch(U.GENCODE_V40_LIFT37, str(GTF) + ".gz")
    subprocess.check_call(["gunzip", "-f", str(GTF) + ".gz"])

# rMATS runs with this kernel's Python (or, on Colab, from its own environment: RMATS_CMD, set by colab_setup)
RMATS = os.environ["RMATS_CMD"].split() if "RMATS_CMD" in os.environ else [sys.executable, os.path.realpath(shutil.which("rmats.py"))]
version = subprocess.run(RMATS + ["--version"], capture_output=True, text=True).stdout.strip()
settings = "annotation {}\\nlibType {}\\nrMATS {}".format(GTF.name, U.RMATS_SETTINGS["lib_type"], version)
print(settings)

manifest = pd.read_csv(BAM_DIR / ("manifest_qki_lite.tsv" if U.LITE else "manifest_qki.tsv"), sep="\\t")
se_tables = {}
for key in U.KNOCKDOWNS:
    print(key, "running rMATS..." if not (DIRS["rmats_run"] / key / "SE.MATS.JC.txt").exists() else "", flush=True)
    se_tables[key], status = U.run_rmats(key, manifest, GTF, DIRS["rmats_run"], THREADS, RMATS, settings)
    print(key, status, se_tables[key])"""),
    ("md", """In lite mode, the next cell compares the PICALM rows of this small run with the full run. Junction reads, PSI and ΔPSI match: rMATS counts each exon's reads on their own. The FDR differs, because it corrects for the number of exons tested, here a handful instead of tens of thousands. The rest of the workflow uses the full run's table."""),
    ("code", """if U.LITE:
    rows = []
    for key in U.KNOCKDOWNS:
        for exon, info in U.EXONS.items():
            for source, path in (("lite run (PICALM only)", se_tables[key]), ("full run", DIRS["rmats"] / key / "SE.MATS.JC.txt")):
                st = U.exon_stats(U.exon_row(path, exon))
                rows.append({"knockdown": key, "exon": info["name"], "source": source, "inclusion reads (KD)": st["IJC_kd"],
                             "skipping reads (KD)": st["SJC_kd"], "PSI KD": st["psi_kd"], "dPSI": st["dpsi"],
                             "P": st["pvalue"], "FDR": st["fdr"]})
    display(pd.DataFrame(rows))
    se_tables = {key: DIRS["rmats"] / key / "SE.MATS.JC.txt" for key in U.KNOCKDOWNS}"""),
    ("md", """## Three groups of exons

| Group | Rule |
| --- | --- |
| included upon knockdown | FDR < 0.1, P < 0.05, ΔPSI > 0.05 |
| excluded upon knockdown | FDR < 0.1, P < 0.05, ΔPSI < −0.05 |
| background | P > 0.5, PSI between 0.05 and 0.95 in both conditions |

An exon included more after knockdown is one QKI normally represses; one excluded more is one QKI normally promotes. Background exons are alternative but did not change.

Exons of the genes in `U.EXCLUDED_GENES` (one Ensembl gene ID) are left out of all three groups.

rMATS can report one exon several times, for example with two different 3′ splice sites. Each group keeps one event per exon by rbp-maps' rule: the form with the most inclusion reads. When two forms tie, the one listed first in the rMATS table is kept, so every computer keeps the same events. Then all three exons of each event (the cassette exon and its two neighbours) are converted from hg19 to GRCh38 with UCSC's chain file, and events that do not convert cleanly are dropped."""),
    ("code", """rows = []
for key in U.KNOCKDOWNS:
    counts, lists = U.build_event_lists(key, se_tables[key], EVENT_DIR, GENOME_CACHE)
    for set_name, n in counts.items():
        rows.append({"knockdown": key, "group": set_name, **n})
    for exon, info in U.EXONS.items():
        print(key, info["name"], "->", U.which_set(lists, exon))
event_counts = pd.DataFrame(rows)
event_counts.to_csv(RESULTS / "05_event_counts.tsv", sep="\\t", index=False)
event_counts"""),
    ("md", """Expected: in both knockdowns the 150-nt exon is in the included group and the 24-nt exon in the excluded group. For HepG2, about 1,500 included, 900 excluded and 2,300 background exons remain."""),
]

# ------------------------------------------------------------------ 06 dPSI per exon
NB["06_dpsi_per_exon"] = [
    ("md", """# 06 · ΔPSI per exon

Notebook 05 tested every exon. This notebook pulls out the two PICALM exons: junction reads, PSI per replicate, ΔPSI and FDR, together with the QKI knockdown level from notebook 04. These numbers make panel d of the figure.

| Exon | GENCODE ID | hg19 (rMATS) | GRCh38 |
| --- | --- | --- | --- |
| 150-nt exon | ENSE00000989276.1 | chr11:85,701,292-85,701,442 | chr11:85,990,249-85,990,399 |
| 24-nt exon | ENSE00001416869.1 | chr11:85,689,112-85,689,136 | chr11:85,978,069-85,978,093 |

Both are on the minus strand; starts are 0-based, as in rMATS."""),
    ("code", SETUP + '''
RMATS_DIR = U.data_dirs(DATA)["rmats"]       # in lite mode: the full run's tables in lite_data/'''),
    ("code", """level = pd.read_csv(RESULTS / "04_knockdown_level.tsv", sep="\\t").set_index("knockdown")
knockdown, rows = {}, []
for key in U.KNOCKDOWNS:
    knockdown[key] = {"log2fc": float(level.loc[key, "QKI log2FC"])}
    for exon, info in U.EXONS.items():
        row = U.exon_row(RMATS_DIR / key / "SE.MATS.JC.txt", exon)
        st = U.exon_stats(row)
        knockdown[key][exon] = st
        rows.append({"knockdown": key, "exon": info["name"], "PSI control": st["psi_ctl"], "PSI knockdown": st["psi_kd"],
                     "dPSI": st["dpsi"], "FDR": st["fdr"], "inclusion reads (ctl | KD)": (st["IJC_ctl"], st["IJC_kd"]),
                     "skipping reads (ctl | KD)": (st["SJC_ctl"], st["SJC_kd"]), "note": row["note"]})
(RESULTS / "06_picalm_knockdown.json").write_text(json.dumps(knockdown, indent=1))
pd.DataFrame(rows)"""),
    ("code", """fig = U.plot_psi(knockdown)"""),
    ("md", """QKI knockdown raises the 150-nt exon's PSI (+0.22 in K562, +0.47 in HepG2) and lowers the 24-nt exon's (−0.10 and −0.55). Both changes are larger in HepG2, where less QKI is left. FDR is below 1e-10 in every case."""),
]

# ------------------------------------------------------------------ 07 integrate
NB["07_integrate"] = [
    ("md", """# 07 · Integrate

This notebook puts the pieces together and draws the figure. It answers three questions:

1. Which transcripts carry each exon, and where does QKI bind next to it? (panels b and c, with notebook 02's windows)
2. How does QKI knockdown change each exon? (panel d, from notebook 06)
3. Do the two PICALM exons follow QKI's rule across the genome? (panel e, splice maps)"""),
    ("code", SETUP + '''
import numpy as np

GENOME_CACHE = DATA / "downloads" / "genome"
DIRS = U.data_dirs(DATA)
EVENT_DIR, MAP_DIR, SIGNAL_DIR = DIRS["events"], DIRS["maps"], DIRS["signal"]
FIGURE = ROOT / "figures" / "workflow" / "CS3_PICALM_QKI_workflow"
import os'''),
    ("md", """## 1. Isoforms and binding

GENCODE v40 lists PICALM's transcripts. Among the protein-coding ones, how many include or skip each exon?"""),
    ("code", """transcripts = U.gencode_transcripts(U.GENE, GENOME_CACHE)
isoforms = U.isoform_table(transcripts)
coding = isoforms[isoforms["type"] == "protein_coding"]
coding.drop(columns="type")"""),
    ("code", """coding[[c for c in coding.columns if c.endswith(("included", "skipped"))]].sum().rename("coding transcripts")"""),
    ("md", """The MANE transcript (PICALM-202) includes the 150-nt exon and skips the 24-nt exon; no GENCODE transcript has both, though brain cDNA clones with both exist (Parikh et al. 2014).

QKI binds the sequence ACUAA, best as ACUAAY with a second UAAY a few nucleotides away (Galarneau & Richard 2005). The next cell counts these motifs in the 100 nt of intron next to each exon, on the RNA strand."""),
    ("code", """import re

chrom, strand = U.LOCUS[:2]
e13, e24 = U.EXONS["exon13"]["grch38"][2:], U.EXONS["exon24nt"]["grch38"][2:]
regions = {"before the 150-nt exon": (e13[1], e13[1] + 100), "after the 24-nt exon": (e24[0] - 100, e24[0])}
for label, (start, end) in regions.items():
    rna = U.rna_sequence(chrom, start, end, strand, GENOME_CACHE)
    print(f"{label:24s} ACUAAY: {len(re.findall('ACUAA[CU]', rna))}  UAAY: {len(re.findall('UAA[CU]', rna))}")
    print("   ", rna)"""),
    ("md", """## 2. Splice maps

A splice map averages eCLIP signal (IP minus input) across many exons, lined up on their splice sites. If QKI represses exons by binding before them, the exons included more after knockdown should have extra signal just before their 3′ splice site. If QKI promotes exons by binding after them, the exons excluded more should have extra signal just after their 5′ splice site. The gray band shows where 99% of random sets of background exons fall.

Two maps use the HepG2 knockdown, the stronger one: one with the ENCORE2 eCLIP (made in K562) and one with the ENCODE3 eCLIP from HepG2 itself. Each replicate takes about 10 minutes with rbp-maps' permutation test. Finished maps are read back from disk, so this cell needs the BAMs from notebook 03 only for maps that are missing."""),
    ("code", """REPS = os.environ.get("REPS", "rep1,rep2").split(",")     # set REPS=rep1 to compute replicate 1 only
maps, tracks = {}, {}
for subject, kd_key in U.MAP_PAIRS:
    lists = U.event_lists(kd_key, EVENT_DIR)
    for rep in REPS:
        if U.map_is_finished(MAP_DIR, subject, kd_key, rep, lists):
            print(subject, kd_key, rep, "-> finished earlier, reading it back")
        else:
            print(subject, kd_key, rep, "-> computing (about 10 min)", flush=True)
            bams = pd.read_csv(RESULTS / "03_eclip_bams.tsv", sep="\\t").set_index(["subject", "sample"])
            ip_sample, input_sample = U.ECLIP[subject]["reps"][rep]
            for sample in (ip_sample, input_sample):
                if (subject, sample) not in tracks:
                    row = bams.loc[(subject, sample)]
                    library = row[["mapped_F260", "mapped_records"]].to_dict() if U.LITE else None
                    tracks[(subject, sample)] = U.signal_tracks(row["bam"], lists[U.SETS[2]], SIGNAL_DIR, library=library)
            U.run_splice_map(tracks[(subject, ip_sample)], tracks[(subject, input_sample)], lists,
                             U.map_png(MAP_DIR, subject, kd_key, rep))
        maps[(subject, rep)] = U.load_splice_map(MAP_DIR, subject, kd_key, rep, lists)
display(Image(str(U.map_png(MAP_DIR, U.MAP_PAIRS[0][0], U.MAP_PAIRS[0][1], "rep1")), width=800))"""),
    ("md", """The plot above is rbp-maps' own. Below, each PICALM exon's own row of the same data: where its eCLIP signal peaks, and how many positions of the genome-wide lines lie above the gray band."""),
    ("code", """REGION = {"exon13": "3ss", "exon24nt": "5ss"}     # QKI binds before the 150-nt exon and after the 24-nt exon
genome_wide, picalm = [], []
for (subject, rep), m in maps.items():
    kd_key = dict(U.MAP_PAIRS)[subject]
    genome_wide.append({"eCLIP": subject, "rep": rep,
                        "included exons: positions above band, 100 nt before 3'ss": U.positions_above_band(
                            m["means"][U.SETS[0]], m["band"][1], "3ss"),
                        "excluded exons: positions above band, 100 nt after 5'ss": U.positions_above_band(
                            m["means"][U.SETS[1]], m["band"][1], "5ss")})
    m["traces"] = {}
    for exon, info in U.EXONS.items():
        trace, set_name = U.exon_trace(m["matrices"], exon)
        m["traces"][exon] = trace
        np.savetxt(MAP_DIR / f"{kd_key}_KD" / f"PICALM_{exon}_{subject}_{rep}.txt", trace)
        picalm.append({"eCLIP": subject, "rep": rep, "exon": info["name"], "group": set_name,
                       "peak (nt)": U.peak(trace, REGION[exon]), "max (x1e-4)": round(float(np.nanmax(trace)) * 1e4, 1)})
genome_wide, picalm = pd.DataFrame(genome_wide), pd.DataFrame(picalm)
genome_wide.to_csv(RESULTS / "07_splice_maps_genome_wide.tsv", sep="\\t", index=False)
picalm.to_csv(RESULTS / "07_splice_maps_picalm.tsv", sep="\\t", index=False)
display(genome_wide)
picalm"""),
    ("md", """Expected for replicate 1. With the ENCORE2 eCLIP, the included exons are above the band at 85 of the 100 positions before the 3′ splice site and the excluded exons at 86 of 100 after the 5′ splice site; with the HepG2 eCLIP, 54 and 86. PICALM's own signal peaks 33 nt (ENCORE2) or 11 nt (HepG2) before the 150-nt exon, and 60 or 68 nt after the 24-nt exon. Both exons sit where the genome-wide rule says QKI binds."""),
    ("md", """## 3. The figure

`draw_figure()` (in `figure.py`) takes the windows from notebook 02, the transcripts above, the numbers from notebook 06 and the replicate-1 maps."""),
    ("code", """from figure import draw_figure

windows = pd.read_csv(RESULTS / "02_qki_windows_PICALM.tsv", sep="\\t")
locus = {"exons": {exon: list(info["grch38"][2:]) for exon, info in U.EXONS.items()}, "transcripts": transcripts,
         "windows": windows[["subject", "start", "end", "log2_enrichment", "q_value"]].values.tolist()}
knockdown = json.loads((RESULTS / "06_picalm_knockdown.json").read_text())
labels = {"QKI_K562_5007": "QKI ENCORE2\\n(K562) eCLIP ×\\nHepG2 KD", "QKI_HepG2_ENCSR570WLM": "QKI ENCODE3\\n(HepG2) eCLIP ×\\nHepG2 KD"}
rows_e = [{"label": labels[subject], "means": maps[(subject, "rep1")]["means"], "band": maps[(subject, "rep1")]["band"],
           "traces": maps[(subject, "rep1")]["traces"]} for subject, _ in U.MAP_PAIRS]
png = draw_figure(locus, knockdown, rows_e, FIGURE)
print("saved", png, "and .pdf")
display(Image(png, width=900))"""),
]

NB["00_colab_run_all"] = [
    ("md", """# 00 · Run all notebooks on Google Colab

[Open in Colab](https://colab.research.google.com/github/yeolab/encode-rbp-case-study/blob/main/notebooks/workflow/00_colab_run_all.ipynb), then choose Runtime > Run all. Nothing needs to be uploaded: the first code cell copies the repository, with its small data, onto the Colab machine.

This notebook runs notebooks 01 to 07 one after another on this Colab machine and saves each finished notebook, with all its tables and plots, to `My Drive/CS3_workflow/` (`executed/`), together with the figure (`figures/workflow/`).

With `LITE = True` (the default, in the next cell) it uses the small data in `lite_data/`: the BAM files cut down to the reads these analyses use. That takes about an hour and fits the free Colab plan. Output goes to `colab_output_lite/`.

With `LITE = False` it downloads everything (about 50 GB) and reruns rMATS genome-wide, which takes 5 to 8 hours. The large files then live on this Colab machine, which is why all seven notebooks run here, in one session. Output goes to `colab_output/`. Keep the browser tab open, or use a paid Colab plan that keeps running in the background."""),
    ("md", """Options: `REP1_ONLY = True` computes only replicate 1 of the splice maps (the one the figure shows), which saves about 40 minutes; it is the default in lite mode. Set `START` to a later notebook to skip the ones already done."""),
    ("code", """import os
import subprocess
import time
from pathlib import Path

REP1_ONLY = os.environ.get("CS3_LITE") == "1"     # lite mode: replicate 1 only, to stay near one hour
START = "01"
NOTEBOOKS = ["01_gene_search", "02_experiment_detail", "03_download", "04_shrna_rnaseq",
             "05_differential_splicing", "06_dpsi_per_exon", "07_integrate"]

if REP1_ONLY:
    os.environ["REPS"] = "rep1"
executed = Path(os.environ.get("CASE_STUDY_ROOT", "../..")).resolve() / "executed"
executed.mkdir(parents=True, exist_ok=True)
for name in NOTEBOOKS:
    if name[:2] < START:
        continue
    started = time.time()
    print(time.strftime("%H:%M"), "running", name, flush=True)
    run = subprocess.run(["jupyter", "nbconvert", "--to", "notebook", "--execute", name + ".ipynb",
                          "--output-dir", str(executed), "--ExecutePreprocessor.timeout=-1",
                          "--ExecutePreprocessor.kernel_name=python3"], capture_output=True, text=True)
    if run.returncode != 0:
        print(run.stderr[-3000:])
        raise RuntimeError(name + " stopped with an error (above). Fix it, then set START = '" + name[:2] + "' and run this cell again.")
    print("      finished in {:.0f} min: {}".format((time.time() - started) / 60, executed / (name + ".ipynb")), flush=True)
print("Figure:", Path(os.environ.get("CASE_STUDY_ROOT", "../..")).resolve() / "figures" / "workflow" / "CS3_PICALM_QKI_workflow.png")"""),
    ("code", """from IPython.display import Image

Image(str(Path(os.environ.get("CASE_STUDY_ROOT", "../..")) / "figures" / "workflow" / "CS3_PICALM_QKI_workflow.png"), width=900)"""),
]

TITLES = {
    "00_colab_run_all": "00_colab_run_all.ipynb",
    "01_gene_search": "01_gene_search.ipynb", "02_experiment_detail": "02_experiment_detail.ipynb",
    "03_download": "03_download.ipynb", "04_shrna_rnaseq": "04_shrna_rnaseq.ipynb",
    "05_differential_splicing": "05_differential_splicing.ipynb", "06_dpsi_per_exon": "06_dpsi_per_exon.ipynb",
    "07_integrate": "07_integrate.ipynb",
}


def build(out_dir):
    """Writes notebooks 00-07 into `out_dir`, each with the Colab setup cell and its screenshots attached."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for key, cells in NB.items():
        extra = ", rmats=True" if key in ("05_differential_splicing", "00_colab_run_all") else ""
        cells = [cells[0], ("md", COLAB_MD), ("code", COLAB.format(extra=extra))] + cells[1:]
        nb = nbformat.v4.new_notebook()
        nb.metadata = {"kernelspec": {"name": "python3", "display_name": "Python 3 (encode-rbp-case-study)", "language": "python"},
                       "language_info": {"name": "python"}}
        for cell in cells:
            if cell[0] == "shot":
                _, name, alt = cell
                mime = "image/jpeg" if name.endswith(".jpg") else "image/png"
                data = base64.b64encode((SHOTS / name).read_bytes()).decode()
                c = nbformat.v4.new_markdown_cell(f"![{alt}](attachment:{name})")
                c["attachments"] = {name: {mime: data}}
                nb.cells.append(c)
            elif cell[0] == "md":
                nb.cells.append(nbformat.v4.new_markdown_cell(cell[1].strip()))
            else:
                nb.cells.append(nbformat.v4.new_code_cell(cell[1].strip()))
        nbformat.validate(nb)
        for i, c in enumerate(nb.cells):
            if c.cell_type == "code":
                compile(c.source, f"{key}[{i}]", "exec")
        nbformat.write(nb, out_dir / TITLES[key])
        print("wrote", out_dir / TITLES[key], len(nb.cells), "cells")


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else HERE / "out")
