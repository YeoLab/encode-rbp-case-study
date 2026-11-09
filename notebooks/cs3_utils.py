"""Shared settings and helpers for the Case Study 3 notebooks (QKI and two PICALM cassette exons).

Data sources, all public:
  RBP-ARK data API and bulk downloads (Skipper 2.0.0): {RBPARK}/api/downloads/experiments.tsv, windows.tsv.gz
  ENCODE portal: knockdown rMATS ("differential splicing quantifications", hg19) and DESeq2 tables
  UCSC: hg19ToHg38 liftover chain; API for sequence
  GENCODE v40 (GRCh38) for transcripts; GENCODE v40 GRCh37 mapping for rMATS (notebook 02)
"""
import gzip
import io
import json
import os
import tarfile
import time
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
import requests

RBPARK = os.environ.get("RBPARK_API", "https://rbp-ark.com")   # RBP-ARK production site and public data API
ENCODE = "https://www.encodeproject.org"
UCSC_API = "https://api.genome.ucsc.edu"
CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
ALIASES = {"TDP43": "TARDBP", "FMRP": "FMR1", "PRK2": "PKN2", "PARP": "PARP1", "RECQ1": "RECQL", "THOC4": "ALYREF"}
SE_COLS = ["geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE",
           "downstreamES", "downstreamEE"]
HEADERS = {"Accept": "application/json", "User-Agent": "rbpark-case-study-notebooks"}


def fetch(url, dest, params=None, retries=3):
    """Download once to dest (reused on later runs)."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        r = requests.get(url, params=params, stream=True, timeout=600, headers=HEADERS)
        if r.status_code == 503 and attempt < retries - 1:      # RBP-ARK builds windows.tsv.gz on demand
            time.sleep(int(r.headers.get("Retry-After", 60)))
            continue
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as f:
            for block in r.iter_content(1 << 22):
                f.write(block)
        tmp.rename(dest)
        return dest
    raise RuntimeError("could not download " + url)


def encode_json(path, **params):
    """Returns the ENCODE portal object at `path` (for example /files/ENCFF096LJD/) as JSON; `params` are extra query
    parameters.
    """
    params.setdefault("format", "json")
    r = requests.get(ENCODE + path, params=params, headers=HEADERS, timeout=120)
    r.raise_for_status()
    return r.json()


# ---------------- RBP-ARK ----------------
def rbpark_experiments(cache):
    """Returns RBP-ARK's experiment catalog (experiments.tsv) as a table, downloaded once into `cache`."""
    return pd.read_csv(fetch(RBPARK + "/api/downloads/experiments.tsv", Path(cache) / "experiments.tsv"), sep="\t")


def rbpark_windows(cache, rbps=None, chunksize=1_000_000):
    """RBP-ARK's table of reproducible enriched windows (one row per window x experiment), optionally only some RBPs."""
    path = fetch(RBPARK + "/api/downloads/windows.tsv.gz", Path(cache) / "windows.tsv.gz")
    parts = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=chunksize,
                             dtype={"chr": str, "strand": str, "gene_name": str, "gene_id": str}):
        if rbps is not None:
            chunk = chunk[chunk["rbp"].isin(rbps)]
        parts.append(chunk)
    return pd.concat(parts, ignore_index=True)


# ---------------- ENCODE knockdown ----------------
def knockdown_files(output_type):
    """All released ENCODE knockdown files of one output type, with target and cell line."""
    g = encode_json("/search/", type="File", output_type=output_type, status="released", limit="all",
                    field=["accession", "dataset", "file_size", "target.label", "biosample_ontology.term_name"])["@graph"]
    return pd.DataFrame([{"accession": f["accession"], "experiment": f["dataset"].split("/")[2],
                          "target": f.get("target", {}).get("label"), "cell": f["biosample_ontology"]["term_name"],
                          "size": f.get("file_size", 0)} for f in g])


def rmats_se(accession, cache):
    """SE.MATS.JunctionCountOnly.txt from the 'MATS_output' tar (not MATS_Norm_output), or None."""
    path = fetch(f"{ENCODE}/files/{accession}/@@download/{accession}.tar.gz", Path(cache) / f"{accession}.tar.gz")
    with tarfile.open(path) as tar:
        for m in tar.getmembers():
            if m.name.endswith("/MATS_output/SE.MATS.JunctionCountOnly.txt"):
                df = pd.read_csv(tar.extractfile(m), sep="\t")
                df["geneSymbol"] = df["geneSymbol"].astype(str).str.strip('"')
                return df
    return None


def one_event_per_exon(rows, col):
    """Returns the rMATS IDs of the events of `rows` (split rMATS lines; `col` maps column names to indices) kept
    as one event per exon, by the rule of rbp-maps' subset_rmats_junctioncountonly.py: events whose spans from the
    upstream exon's end to the downstream exon's start overlap or touch (on either strand) are merged, and the one
    with the highest mean inclusion junction count is kept. rbp-maps leaves ties to bedtools sort, whose order
    differs between macOS and Linux; here a tie goes to the event listed first, so every computer keeps the same events.
    """
    def mean_ijc(r):
        counts = [int(x) for c in ("IJC_SAMPLE_1", "IJC_SAMPLE_2") for x in r[col[c]].split(",")]
        return sum(counts) / len(counts)
    spans = sorted((r[col["chr"]], int(r[col["upstreamEE"]]), int(r[col["downstreamES"]]), i) for i, r in enumerate(rows))
    keep, group, chrom, end = [], [], None, None
    for c, s, e, i in spans + [(None, 0, 0, -1)]:
        if group and (c != chrom or s > end):
            best = max(group, key=lambda j: (mean_ijc(rows[j]), -j))
            keep.append(rows[best][0])
            group = []
        if c is None:
            break
        if not group:
            chrom, end = c, e
        group.append(i)
        end = max(end, e)
    return keep


def gene_ids(table):
    """Returns the Ensembl gene ID of each row of rMATS table `table`, without quotes or version."""
    return table["GeneID"].astype(str).str.strip('"').str.split(".").str[0]


def se_for(target, cell, files, cache):
    """The SE table for one knockdown (tries the target/cell's rMATS tars, smallest first)."""
    for acc in files[(files.target == target) & (files.cell == cell)].sort_values("size")["accession"]:
        df = rmats_se(acc, cache)
        if df is not None:
            return acc, df[~gene_ids(df).isin(EXCLUDED_GENES)].reset_index(drop=True)
    return None, None


def add_se_stats(df):
    """Junction reads, mean PSI per condition. rMATS sample 1 = knockdown, sample 2 = control."""
    split_sum = lambda s: s.astype(str).str.split(",").apply(lambda v: sum(float(x) for x in v if x not in ("", "NA")))
    split_mean = lambda s: s.astype(str).str.split(",").apply(
        lambda v: np.nanmean([float(x) for x in v if x not in ("", "NA")]) if any(x not in ("", "NA") for x in v) else np.nan)
    df = df.copy()
    df["junction_reads"] = sum(split_sum(df[c]) for c in ["IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2"])
    df["psi_kd"], df["psi_ctl"] = split_mean(df["IncLevel1"]), split_mean(df["IncLevel2"])
    df["dpsi"] = df["IncLevelDifference"]
    return df


def event_sets(df, min_dpsi=0.05):
    """ENCODE thresholds: included, excluded upon knockdown, and native cassette background."""
    df = add_se_stats(df)
    ok = df["junction_reads"] >= 30
    sig = ok & (df["FDR"] < 0.05) & (df["PValue"] < 0.05)
    return {"inc": df[sig & (df["dpsi"] >= min_dpsi)], "exc": df[sig & (df["dpsi"] <= -min_dpsi)],
            "bg": df[ok & (df["FDR"] >= 0.5) & df["psi_ctl"].between(0.05, 0.95)]}


def deseq2_log2fc(experiment, gene_id, cache):
    """Knockdown level in each of an experiment's DESeq2 tables, smallest file first."""
    g = encode_json("/search/", type="File", dataset=f"/experiments/{experiment}/",
                    output_type="differential expression quantifications", status="released",
                    field=["accession", "file_size"])["@graph"]
    out = {}
    for f in sorted(g, key=lambda x: x["file_size"]):
        acc = f["accession"]
        try:
            t = pd.read_csv(fetch(f"{ENCODE}/files/{acc}/@@download/{acc}.tsv", Path(cache) / f"{acc}.tsv"),
                            sep="\t", index_col=0)
        except Exception as e:
            print("skipped", acc, e)
            continue
        if "log2FoldChange" in t.columns:
            hit = t[t.index.astype(str).str.startswith(gene_id)]
            if len(hit):
                out[acc] = float(hit["log2FoldChange"].iloc[0])
    return out


# ---------------- liftover ----------------
_LIFT = None


def lifter(cache):
    """Returns a pyliftover hg19 -> GRCh38 converter, built once from UCSC's chain file (downloaded into `cache`)."""
    global _LIFT
    if _LIFT is None:
        from pyliftover import LiftOver
        _LIFT = LiftOver(str(fetch(CHAIN_URL, Path(cache) / "hg19ToHg38.over.chain.gz")))
    return _LIFT


def lift_events(df, cache, cols=("exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE")):
    """Lift each event's coordinates to GRCh38; keep events whose coordinates all lift to the same chromosome
    with the exon length unchanged."""
    lo = lifter(cache)
    rows = []
    for _, r in df.iterrows():
        new = {}
        for c in cols:
            m = lo.convert_coordinate(r["chr"], int(r[c]))
            if not m or m[0][0] != r["chr"]:
                break
            new[c] = m[0][1]
        else:
            if new["exonEnd"] - new["exonStart_0base"] == int(r["exonEnd"]) - int(r["exonStart_0base"]):
                rows.append({**r.to_dict(), **{c + "_hg38": v for c, v in new.items()}})
    return pd.DataFrame(rows)


# ---------------- windows near exons ----------------
class WindowIndex:
    """Bound windows (q <= qmax, log2 >= min_l2) per experiment, chromosome and strand, for overlap queries."""

    def __init__(self, windows, qmax=0.05, min_l2=2.0):
        """Indexes the windows in `windows` with q <= `qmax` and log2 enrichment >= `min_l2` by experiment, chromosome
        and strand.
        """
        w = windows[(windows["q_value"] <= qmax) & (windows["log2_enrichment"] >= min_l2)]
        self.idx = {k: (g["start"].to_numpy(), g["end"].to_numpy())
                    for k, g in w.sort_values("start").groupby(["subject", "chr", "strand"])}

    def any(self, subject, chrom, strand, s, e):
        """Returns True if experiment `subject` has a bound window overlapping `chrom`:`s`-`e` on `strand`."""
        hit = self.idx.get((subject, chrom, strand))
        if hit is None:
            return False
        starts, ends = hit
        i = np.searchsorted(starts, e)            # windows starting before e
        lo = max(0, i - 50)
        return bool(np.any(ends[lo:i] > s))


def exon_bins(strand, es, ee, step=50, flank=300):
    """13 bins in transcript orientation: 6 upstream-intron bins (300..0 nt before the 3' splice site), the exon,
    6 downstream-intron bins (0..300 nt after the 5' splice site)."""
    up = [(es - a, es - a + step) if strand == "+" else (ee + a - step, ee + a) for a in range(flank, 0, -step)]
    dn = [(ee + a, ee + a + step) if strand == "+" else (es - a - step, es - a) for a in range(0, flank, step)]
    return up + [(es, ee)] + dn


BIN_NAMES = ["U300", "U250", "U200", "U150", "U100", "U50", "EX", "D50", "D100", "D150", "D200", "D250", "D300"]


def window_map(events, index, subject):
    """Number of exons with a bound window in each bin."""
    counts = np.zeros(len(BIN_NAMES), dtype=int)
    for _, r in events.iterrows():
        for i, (s, e) in enumerate(exon_bins(r["strand"], r["exonStart_0base_hg38"], r["exonEnd_hg38"])):
            counts[i] += index.any(subject, r["chr"], r["strand"], s, e)
    return counts


# ---------------- UCSC ----------------
def ucsc_sequence(chrom, start, end, strand="+"):
    """Returns the GRCh38 sequence of `chrom`:`start`-`end` (0-based, end exclusive) from the UCSC API, as RNA, reverse-
    complemented when `strand` is "-".
    """
    r = requests.get(f"{UCSC_API}/getData/sequence", params={"genome": "hg38", "chrom": chrom, "start": start, "end": end},
                     timeout=60)
    r.raise_for_status()
    s = r.json()["dna"].upper()
    if strand == "-":
        s = s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
    return s.replace("T", "U")


GENCODE_V40 = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz"


def gencode_transcripts(gene, cache):
    """Transcripts of one gene in GENCODE v40 (GRCh38): {transcript_id: {name, type, tags, exons}}, exons 0-based half-open."""
    import re
    path = fetch(GENCODE_V40, Path(cache) / "gencode.v40.annotation.gtf.gz")
    out = {}
    key = 'gene_name "{}";'.format(gene)
    with gzip.open(path, "rt") as f:
        for line in f:
            if key not in line:
                continue
            p = line.rstrip("\n").split("\t")
            if p[2] not in ("transcript", "exon"):
                continue
            attr = p[8]
            tid = re.search(r'transcript_id "([^"]+)"', attr).group(1)
            t = out.setdefault(tid, {"name": re.search(r'transcript_name "([^"]+)"', attr).group(1),
                                     "type": re.search(r'transcript_type "([^"]+)"', attr).group(1),
                                     "tags": re.findall(r'tag "([^"]+)"', attr), "chrom": p[0], "strand": p[6], "exons": []})
            if p[2] == "exon":
                t["exons"].append([int(p[3]) - 1, int(p[4])])
    for t in out.values():
        t["exons"].sort()
    return out


# ---------------- folders ----------------
DATA_DEFAULT = Path(__file__).resolve().parents[1] / "data"      # <case study folder>/data


def folders():
    """Returns (ROOT, DATA). ROOT: the case study folder (notebooks/, figures/, results/), set by CASE_STUDY_ROOT.
    DATA: downloads and other large intermediate files, set by CASE_STUDY_DATA, by default <ROOT>/data."""
    root = Path(os.environ.get("CASE_STUDY_ROOT", "..")).resolve()
    data = Path(os.environ.get("CASE_STUDY_DATA", root / "data"))
    data.mkdir(parents=True, exist_ok=True)
    for d in (root / "figures" / "data", root / "results"):
        d.mkdir(parents=True, exist_ok=True)
    return root, data.resolve()


# ---------------- the four knockdowns (ENCODE shRNA RNA-seq, hg19) ----------------
# rmats: ENCODE's rMATS skipped-exon release. deseq2: the DESeq2 table without batch normalization, used for the
# knockdown level; deseq2_batch: the batch-normalized table of the same comparison, shown for comparison.
KNOCKDOWNS = OrderedDict([
    ("QKI_K562", dict(rbp="QKI", cell="K562", kd_experiment="ENCSR256PLH", control_experiment="ENCSR344XID",
                      rmats="ENCFF096LJD", deseq2="ENCFF834IQO", deseq2_batch="ENCFF557AJB",
                      gene_id="ENSG00000112531")),
    ("QKI_HepG2", dict(rbp="QKI", cell="HepG2", kd_experiment="ENCSR330YOU", control_experiment="ENCSR603TCV",
                       rmats="ENCFF103PLX", deseq2="ENCFF526MSX", deseq2_batch="ENCFF279VLA",
                       gene_id="ENSG00000112531")),
])

# ---------------- the two PICALM exons (minus strand) ----------------
COORDS = ["chr", "strand", "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"]
PICALM_HG19 = OrderedDict([
    ("exon13", dict(chr="chr11", strand="-", exonStart_0base=85701292, exonEnd=85701442, upstreamES=85694908,
                    upstreamEE=85695016, downstreamES=85707868, downstreamEE=85707972)),
    ("exon24nt", dict(chr="chr11", strand="-", exonStart_0base=85689112, exonEnd=85689136, upstreamES=85687665,
                      upstreamEE=85687725, downstreamES=85692171, downstreamEE=85692271)),
])
PICALM_GRCH38 = {"exon13": ("chr11", "-", 85990249, 85990399), "exon24nt": ("chr11", "-", 85978069, 85978093)}

# Exons of these genes (Ensembl IDs, without version) are left out of every event list.
EXCLUDED_GENES = {"ENSG00000151923"}
# rbp-maps draws its background subsets with Python's random module; a fixed seed makes the band reproducible.
PERMUTATION_SEED = 0


# ---------------- the three QKI eCLIP datasets on RBP-ARK (Skipper 2.0.0, GRCh38) ----------------
# All three give the windows in figure panels b and c; the two in MAP_PAIRS also give the splice maps in panel e.
ECLIP = OrderedDict([
    ("QKI_K562_5007", dict(rbp="QKI", cell="K562", project="ENCORE2", reps={"rep1": ("IP_1", "IN_1"), "rep2": ("IP_2", "IN_2")})),
    ("QKI_K562_ENCSR366YOG", dict(rbp="QKI", cell="K562", project="ENCODE3", reps={"rep1": ("IP_1", "IN_1"), "rep2": ("IP_2", "IN_1")})),
    ("QKI_HepG2_ENCSR570WLM", dict(rbp="QKI", cell="HepG2", project="ENCODE3", reps={"rep1": ("IP_1", "IN_1"), "rep2": ("IP_2", "IN_1")})),
])

# The splice maps in figure panel e: (eCLIP dataset, knockdown). Both use the HepG2 knockdown, the stronger one.
MAP_PAIRS = [("QKI_K562_5007", "QKI_HepG2"), ("QKI_HepG2_ENCSR570WLM", "QKI_HepG2")]


def rbpark_files(subject, cache):
    """Every file RBP-ARK lists for one dataset (its API record): url, path on the collection, size, sha256."""
    path = fetch(f"{RBPARK}/api/experiments/{subject}/", Path(cache) / f"{subject}.json")
    rows = []

    def walk(x):
        """Collects every file entry (a key/value list with a path) found anywhere in the nested record `x` into `rows`.
        """
        if isinstance(x, list):
            if x and all(isinstance(y, dict) and "key" in y for y in x):
                d = {y["key"]: y.get("value") for y in x}
                if d.get("path"):
                    rows.append(d)
            else:
                for y in x:
                    walk(y)
        elif isinstance(x, dict):
            for y in x.values():
                walk(y)

    walk(json.loads(Path(path).read_text())["content"]["filelist"])
    df = pd.DataFrame(rows).rename(columns={"path": "url", "filesize": "size"})
    df["size"] = df["size"].astype("int64")
    df["collection_path"] = df["url"].str.replace(r"^https://[^/]+", "", regex=True)
    return df[["url", "collection_path", "size", "sha256"]].drop_duplicates("url")
