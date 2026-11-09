"""Functions for the Case Study 3 workflow notebooks (QKI and two PICALM exons).

One notebook per box of figure panel a:

    RBP-ARK (eCLIP)                 01 gene search, 02 experiment detail, 03 download
    ENCODE (shRNA knockdown)        04 shRNA RNA-seq, 05 differential splicing, 06 dPSI per exon
    both                            07 integrate (rbp-maps and the figure)

Every download is cached: a file that is already on disk is not fetched again, so the notebooks can be rerun.
Coordinates follow rMATS and RBP-ARK: 0-based start, half-open end.
"""
import gzip
import hashlib
import json
import os
import random
import shutil
import subprocess
import tarfile
import time
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd
import requests

os.environ.setdefault("TQDM_DISABLE", "1")      # rbp-maps' permutation progress bars would fill the notebook

# ---------------------------------------------------------------- where things are
RBPARK = os.environ.get("RBPARK_API", "https://rbp-ark.com")
ENCODE = "https://www.encodeproject.org"
UCSC_API = "https://api.genome.ucsc.edu"
CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
GENCODE_V40 = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz"
GENCODE_V40_LIFT37 = ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh37_mapping/"
                      "gencode.v40lift37.annotation.gtf.gz")
DATA_DEFAULT = Path(__file__).resolve().parents[2] / "data"      # <case study folder>/data
HEADERS = {"Accept": "application/json", "User-Agent": "rbpark-case-study-notebooks"}

# ---------------------------------------------------------------- the example
GENE = "PICALM"
RBP = "QKI"

# The three QKI eCLIP experiments on RBP-ARK (Skipper 2.0.0, GRCh38). reps: (IP sample, size-matched input sample)
ECLIP = OrderedDict([
    ("QKI_K562_5007", dict(cell="K562", project="ENCORE2", reps={"rep1": ("IP_1", "IN_1"), "rep2": ("IP_2", "IN_2")})),
    ("QKI_K562_ENCSR366YOG", dict(cell="K562", project="ENCODE3", reps={"rep1": ("IP_1", "IN_1"), "rep2": ("IP_2", "IN_1")})),
    ("QKI_HepG2_ENCSR570WLM", dict(cell="HepG2", project="ENCODE3", reps={"rep1": ("IP_1", "IN_1"), "rep2": ("IP_2", "IN_1")})),
])

# The two ENCODE QKI shRNA knockdowns (hg19). rmats: ENCODE's rMATS release; deseq2: DESeq2 table without batch
# normalization (used for the knockdown level); deseq2_batch: the batch-normalized table of the same comparison.
KNOCKDOWNS = OrderedDict([
    ("QKI_K562", dict(cell="K562", kd_experiment="ENCSR256PLH", control_experiment="ENCSR344XID",
                      rmats="ENCFF096LJD", deseq2="ENCFF834IQO", deseq2_batch="ENCFF557AJB", gene_id="ENSG00000112531")),
    ("QKI_HepG2", dict(cell="HepG2", kd_experiment="ENCSR330YOU", control_experiment="ENCSR603TCV",
                       rmats="ENCFF103PLX", deseq2="ENCFF526MSX", deseq2_batch="ENCFF279VLA", gene_id="ENSG00000112531")),
])

# The two PICALM cassette exons (minus strand). hg19 rows are rMATS events (exon plus its two flanking exons).
EXONS = OrderedDict([
    ("exon13", dict(name="150-nt exon", gencode_id="ENSE00000989276.1", grch38=("chr11", "-", 85990249, 85990399),
                    hg19=dict(chr="chr11", strand="-", exonStart_0base=85701292, exonEnd=85701442, upstreamES=85694908,
                              upstreamEE=85695016, downstreamES=85707868, downstreamEE=85707972))),
    ("exon24nt", dict(name="24-nt exon", gencode_id="ENSE00001416869.1", grch38=("chr11", "-", 85978069, 85978093),
                      hg19=dict(chr="chr11", strand="-", exonStart_0base=85689112, exonEnd=85689136, upstreamES=85687665,
                                upstreamEE=85687725, downstreamES=85692171, downstreamEE=85692271))),
])
LOCUS = ("chr11", "-", 85975000, 85995000)          # GRCh38 region holding both exons

# Splice maps drawn in the figure: (eCLIP experiment, knockdown)
MAP_PAIRS = [("QKI_K562_5007", "QKI_HepG2"), ("QKI_HepG2_ENCSR570WLM", "QKI_HepG2")]
SETS = ["included-upon-knockdown", "excluded-upon-knockdown", "native-cassette-background"]
# Exons of these genes (Ensembl IDs, without version) are left out of every event list.
EXCLUDED_GENES = {"ENSG00000151923"}
# rbp-maps draws its background subsets with Python's random module; a fixed seed makes the band reproducible.
PERMUTATION_SEED = 0

# Lite mode (free Google Colab): small BAMs that keep only the reads these analyses use, from lite_data/
# (made once by make_lite_data.py from the full data). Turned on by colab_setup.setup(..., lite=True).
LITE = os.environ.get("CS3_LITE") == "1"
LITE_DIR = Path(os.environ.get("CS3_LITE_DIR", Path(__file__).resolve().parent / "lite_data"))


# ---------------------------------------------------------------- folders and downloads
def folders():
    """Return (root, data).

    root: the case study folder (results/ and figures/ go here). Set CASE_STUDY_ROOT, default two levels above
    the notebooks (notebooks/workflow/ -> CASE_STUDY_3/).
    data: downloads and other large files. Set CASE_STUDY_DATA, default <root>/data.
    """
    root = Path(os.environ.get("CASE_STUDY_ROOT", "../..")).resolve()
    data = Path(os.environ.get("CASE_STUDY_DATA", root / "data"))
    data.mkdir(parents=True, exist_ok=True)
    for d in (root / "results" / "workflow", root / "figures" / "workflow"):
        d.mkdir(parents=True, exist_ok=True)
    return root, data.resolve()


def fetch(url, dest, params=None, retries=4):
    """Download url to dest once and return the path. A file already at dest is reused.

    RBP-ARK builds some tables on request and answers 503 with a Retry-After header until they are ready;
    this waits and tries again.
    """
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        r = requests.get(url, params=params, stream=True, timeout=600, headers=HEADERS)
        if r.status_code in (429, 503) and attempt < retries - 1:
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


def checksum(path, algorithm="sha256", chunk=8 << 20):
    """Hex digest of a file (sha256 or md5)."""
    h = hashlib.new(algorithm)
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download_large(url, out, size, digest=None, algorithm="sha256"):
    """Download a large file with curl (resumes a partial file), then check its size and checksum.

    A file that passed the checksum once gets a small '.<algorithm>ok' marker and is not checked again.
    """
    out = Path(out)
    ok = Path(str(out) + "." + algorithm + "ok")
    if ok.exists() and out.exists() and out.stat().st_size == size:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    if not (out.exists() and out.stat().st_size == size):
        subprocess.check_call(["curl", "-L", "-sS", "--retry", "5", "-C", "-", "-o", str(out), url])
    assert out.stat().st_size == size, "{}: size {} != {}".format(out.name, out.stat().st_size, size)
    if digest:
        assert checksum(out, algorithm) == digest, "{}: checksum mismatch; delete it and run again".format(out.name)
        ok.write_text(digest + "\n")
    return out


# ---------------------------------------------------------------- RBP-ARK (notebooks 01-03)
def rbpark_json(path, cache, **params):
    """GET an RBP-ARK API path (for example '/api/genes/PICALM/experiments/') and cache the JSON answer."""
    name = (path.strip("/").replace("/", "_") + ("_" + "_".join(f"{k}-{v}" for k, v in sorted(params.items()))
                                                  if params else "")).replace(",", "+") + ".json"
    return json.loads(fetch(RBPARK + path, Path(cache) / name, params=params or None).read_text())


def gene_search(gene, cache, min_enrichment=3):
    """What rbp-ark.com's gene search shows: every window on a gene, one row per window and experiment.

    Uses GET /api/genes/<gene>/experiments/?min_enrichment=<log2>. Columns: subject (experiment), rbp, chr, start,
    end, strand, feature_type, log2_enrichment, q_value.
    """
    answer = rbpark_json(f"/api/genes/{gene}/experiments/", cache, min_enrichment=min_enrichment)
    rows = [{"subject": e["subject"], "rbp": e["rbp"], **w} for e in answer["experiments"] for w in e["windows"]]
    return pd.DataFrame(rows)


def experiments_table(cache):
    """RBP-ARK's list of every public experiment: subject, rbp, cell_type, status (QC), project, enriched_windows."""
    return pd.read_csv(fetch(RBPARK + "/api/downloads/experiments.tsv", Path(cache) / "experiments.tsv"), sep="\t")


def experiment_record(subject, cache, fields=None):
    """One experiment's API record (GET /api/experiments/<subject>/). fields: comma-separated dotted paths."""
    params = {"fields": fields} if fields else {}
    return rbpark_json(f"/api/experiments/{subject}/", cache, **params)


def _key_value(items):
    """Returns a list of {"key": ..., "value": ...} entries (`items`) as a dict."""
    return {d["key"]: d.get("value") for d in items}


def read_stats(subject, cache):
    """Read counts per sample from an experiment's QC summary (the table at the top of its RBP-ARK page)."""
    stats = experiment_record(subject, cache, "content.summary_data.basic_stats")["content"]["summary_data"]["basic_stats"]
    table = pd.DataFrame({metric: _key_value(values) for metric, values in stats.items()})
    table.index.name = "sample"
    return table


def file_list(subject, cache):
    """Every file RBP-ARK lists for one experiment: url, collection_path, size, sha256.

    This is the same list as the experiment's BDBag (fetch.txt and manifest-sha256.txt).
    """
    record = experiment_record(subject, cache, "content.filelist")
    rows = [_key_value(f) for f in record["content"]["filelist"]]
    df = pd.DataFrame(rows).rename(columns={"path": "url", "filesize": "size"})
    df["size"] = df["size"].astype("int64")
    df["collection_path"] = df["url"].str.replace(r"^https://[^/]+", "", regex=True)
    return df[["url", "collection_path", "size", "sha256"]].drop_duplicates("url").reset_index(drop=True)


def reproducible_windows(subject, cache):
    """An experiment's reproducible enriched windows file (one row per window), downloaded once.

    log2_enrichment is the file's enrichment_l2or_mean and q_value its q_max: the numbers RBP-ARK's gene search shows.
    """
    files = file_list(subject, cache)
    hit = files[files["collection_path"].str.endswith(f"reproducible_enriched_windows/{subject}.reproducible_enriched_windows.tsv.gz")]
    assert len(hit) == 1, f"{subject}: no reproducible windows file in the file list"
    path = fetch(hit.iloc[0]["url"], Path(cache) / "windows" / Path(hit.iloc[0]["collection_path"]).name)
    df = pd.read_csv(path, sep="\t", dtype={"chr": str})
    return df.rename(columns={"enrichment_l2or_mean": "log2_enrichment", "q_max": "q_value"})


def windows_in(df, chrom, strand, start, end):
    """Windows overlapping chrom:start-end on one strand."""
    return df[(df["chr"] == chrom) & (df["strand"] == strand) & (df["end"] > start) & (df["start"] < end)]


def relative_position(start, end, exon):
    """Where a window sits relative to a PICALM exon, in the direction of transcription (nt).

    The window's nearest nucleotide, numbered like intron positions: -1 is the last intron nt before the exon's
    3' splice site, +1 the first intron nt after its 5' splice site, 0 means the window overlaps the exon.
    """
    chrom, strand, es, ee = EXONS[exon]["grch38"]
    if end > es and start < ee:
        return 0
    if strand == "+":
        return -(es - end + 1) if end <= es else start - ee + 1
    return -(start - ee + 1) if start >= ee else es - end + 1


def pick_bam(files, subject, sample):
    """The deduplicated genome BAM of one sample (IP_1, IN_1, ...): the read-2-only file when there is one."""
    pattern = r"/bams/dedup/genome(?:_R2)?/{}_{}\.genome\.Aligned\.sort\.dedup(?:\.R2)?\.bam$".format(subject, sample)
    hits = files[files["collection_path"].str.contains(pattern)]
    assert len(hits), f"{subject}: no deduplicated BAM for {sample}"
    return hits.sort_values("collection_path", key=lambda s: ~s.str.contains("genome_R2")).iloc[0]


def eclip_bams(cache):
    """The eCLIP BAMs the splice maps need (both replicates of every experiment in MAP_PAIRS), with URL, size, sha256."""
    rows = []
    for subject in OrderedDict.fromkeys(name for name, _ in MAP_PAIRS):
        files = file_list(subject, cache)
        for sample in sorted({s for pair in ECLIP[subject]["reps"].values() for s in pair}):
            row = pick_bam(files, subject, sample)
            rows.append({"subject": subject, "sample": sample, **row.to_dict()})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- ENCODE (notebooks 04-06)
def encode_json(path, cache):
    """GET an ENCODE portal object as JSON (path like /files/ENCFF096LJD/) and cache it."""
    name = path.strip("/").replace("/", "_") + ".json"
    return json.loads(fetch(ENCODE + path, Path(cache) / name, params={"format": "json"}).read_text())


def _object_id(x):
    """Returns the ENCODE @id of `x`, which is either an @id string or an embedded object."""
    return x if isinstance(x, str) else x["@id"]


def knockdown_bams(key, cache):
    """The four BAMs ENCODE's rMATS file for one knockdown was computed from (its derived_from field).

    Two shRNA replicates (rMATS sample 1) and two non-target control replicates (sample 2), hg19 TopHat alignments.
    """
    kd = KNOCKDOWNS[key]
    rmats = encode_json(f"/files/{kd['rmats']}/", cache)
    rows = []
    for parent in rmats["derived_from"]:
        bam = encode_json(_object_id(parent), cache)
        experiment = _object_id(bam["dataset"]).strip("/").split("/")[-1]
        rows.append({"knockdown": key, "cell": kd["cell"],
                     "condition": "knockdown" if experiment == kd["kd_experiment"] else "control",
                     "rep": bam["biological_replicates"][0], "experiment": experiment, "accession": bam["accession"],
                     "assembly": bam["assembly"], "read_length": bam.get("mapped_read_length"),
                     "file_size": bam["file_size"], "md5sum": bam["md5sum"], "url": ENCODE + bam["href"]})
    return pd.DataFrame(rows).sort_values(["condition", "rep"], ascending=[False, True]).reset_index(drop=True)


def encode_file(accession, folder, cache):
    """Download one ENCODE file (checked against its MD5) into folder and return the path."""
    meta = encode_json(f"/files/{accession}/", cache)
    out = Path(folder) / Path(meta["href"]).name
    return download_large(ENCODE + meta["href"], out, meta["file_size"], meta["md5sum"], "md5")


def deseq2_log2fc(table_path, gene_id):
    """log2 fold change of one gene (knockdown vs control) in an ENCODE DESeq2 table; NaN if the gene is absent."""
    table = pd.read_csv(table_path, sep="\t", index_col=0)
    hit = [i for i in table.index if str(i) == gene_id or str(i).startswith(gene_id + ".")]
    return float(table.loc[hit[0], "log2FoldChange"]) if hit else float("nan")


def encode_rmats_se(tar_path, out):
    """Pull SE.MATS.JunctionCountOnly.txt (MATS_output) out of an ENCODE rMATS tar."""
    out = Path(out)
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(tar_path) as tar:
            member = next(m for m in tar.getmembers() if m.name.endswith("MATS_output/SE.MATS.JunctionCountOnly.txt"))
            out.write_bytes(tar.extractfile(member).read())
    return out


# ---------------------------------------------------------------- rMATS (notebook 05)
RMATS_SETTINGS = dict(read_length=100, lib_type="fr-firststrand")


def rmats_command(b1, b2, gtf, od, tmp, threads, rmats=("rmats.py",)):
    """The rMATS-turbo command: knockdown BAMs as sample 1, controls as sample 2, so dPSI = knockdown - control."""
    return list(rmats) + ["--b1", str(b1), "--b2", str(b2), "--gtf", str(gtf), "-t", "paired",
                          "--readLength", str(RMATS_SETTINGS["read_length"]), "--variable-read-length",
                          "--libType", RMATS_SETTINGS["lib_type"], "--nthread", str(threads), "--od", str(od),
                          "--tmp", str(tmp)]


def run_rmats(key, manifest, gtf, out_dir, threads, rmats, settings):
    """Run rMATS-turbo for one knockdown, or reuse a finished run made with the same settings.

    out_dir/<key>/SE.MATS.JC.txt is the result. run_settings.txt records the annotation, --libType and rMATS
    version; if they differ from `settings`, the old run is cleared and rMATS runs again.
    """
    sub = manifest[manifest["knockdown"] == key].sort_values("rep")
    kd = sub.loc[sub["condition"] == "knockdown", "bam"].tolist()
    ctl = sub.loc[sub["condition"] == "control", "bam"].tolist()
    assert len(kd) == 2 and len(ctl) == 2, sub
    od, tmp = Path(out_dir) / key, Path(out_dir).parent / "tmp" / key
    od.mkdir(parents=True, exist_ok=True)
    result, stamp = od / "SE.MATS.JC.txt", od / "run_settings.txt"
    if result.exists() and stamp.exists() and stamp.read_text().strip() == settings:
        return result, "reused"
    if result.exists():
        for f in od.iterdir():
            if f.name not in ("b1.txt", "b2.txt"):
                shutil.rmtree(f) if f.is_dir() else f.unlink()
    shutil.rmtree(tmp, ignore_errors=True)         # leftovers from a stopped run make rMATS count BAMs twice
    tmp.mkdir(parents=True, exist_ok=True)
    (od / "b1.txt").write_text(",".join(kd) + "\n")
    (od / "b2.txt").write_text(",".join(ctl) + "\n")
    cmd = rmats_command(od / "b1.txt", od / "b2.txt", gtf, od, tmp, threads, rmats)
    (od / "command.txt").write_text(" ".join(cmd) + "\n")
    with open(od / "rmats.log", "w") as log:
        status = subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT)
    assert status == 0 and result.exists(), "rMATS failed; see " + str(od / "rmats.log")
    stamp.write_text(settings + "\n")
    return result, "ran"


def _floats(text):
    """Returns the numbers in a comma-separated rMATS field `text`, skipping NA."""
    return [float(x) for x in str(text).split(",") if x not in ("NA", "")]


def _mean_psi(text):
    """Returns the mean of the PSI values in rMATS field `text` (NaN if all are NA)."""
    v = _floats(text)
    return float(np.mean(v)) if v else np.nan


def gene_ids(table):
    """Returns the Ensembl gene ID of each row of rMATS table `table`, without quotes or version."""
    return table["GeneID"].astype(str).str.strip('"').str.split(".").str[0]


def classify(table):
    """Boolean masks for the three exon sets of a splice map, from an rMATS SE table (sample 1 = knockdown).

    included-upon-knockdown: FDR < 0.1, P < 0.05, dPSI > 0.05
    excluded-upon-knockdown: FDR < 0.1, P < 0.05, dPSI < -0.05
    native-cassette-background: P > 0.5 and mean PSI strictly between 0.05 and 0.95 in both conditions
    Exons of the genes in EXCLUDED_GENES are in none of the sets.
    """
    psi_kd, psi_ctl = table["IncLevel1"].map(_mean_psi), table["IncLevel2"].map(_mean_psi)
    keep = ~gene_ids(table).isin(EXCLUDED_GENES)
    sig = keep & (table["FDR"] < 0.1) & (table["PValue"] < 0.05)
    return {"included-upon-knockdown": sig & (table["IncLevelDifference"] > 0.05),
            "excluded-upon-knockdown": sig & (table["IncLevelDifference"] < -0.05),
            "native-cassette-background": keep & (table["PValue"] > 0.5)
                                          & psi_kd.between(0.05, 0.95, inclusive="neither")
                                          & psi_ctl.between(0.05, 0.95, inclusive="neither")}


_LIFT = None


def lifter(cache):
    """A pyliftover hg19 -> hg38 converter built from UCSC's chain file (downloaded once)."""
    global _LIFT
    if _LIFT is None:
        from pyliftover import LiftOver
        _LIFT = LiftOver(str(fetch(CHAIN_URL, Path(cache) / "hg19ToHg38.over.chain.gz")))
    return _LIFT


_INTERVALS = [("exonStart_0base", "exonEnd"), ("upstreamES", "upstreamEE"), ("downstreamES", "downstreamEE")]


def lift_event(fields, col, lift):
    """Lift one rMATS line (a list of fields) from hg19 to GRCh38, or return None.

    Kept only if each of its three exons maps whole to one chromosome, on the same strand, with its length
    unchanged (as UCSC liftOver would), and the exons stay in order (upstream < cassette < downstream).
    """
    chrom = fields[col["chr"]]
    out, chroms = list(fields), set()
    for start_name, end_name in _INTERVALS:
        start, end = int(fields[col[start_name]]), int(fields[col[end_name]])
        a, b = lift.convert_coordinate(chrom, start), lift.convert_coordinate(chrom, end - 1)
        if not a or not b or len(a) != 1 or len(b) != 1:
            return None
        a, b = a[0][:3], b[0][:3]
        if a[0] != b[0] or a[2] != "+" or b[2] != "+" or (b[1] + 1) - a[1] != end - start:
            return None
        chroms.add(a[0])
        out[col[start_name]], out[col[end_name]] = str(a[1]), str(b[1] + 1)
    if len(chroms) != 1:
        return None
    out[col["chr"]] = chroms.pop()
    in_order = (int(out[col["upstreamEE"]]) <= int(out[col["exonStart_0base"]])
                and int(out[col["exonEnd"]]) <= int(out[col["downstreamES"]]))
    return out if in_order else None


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


def build_event_lists(key, se_path, out_dir, cache):
    """Write the three GRCh38 event lists rbp-maps reads for one knockdown; return (counts, paths).

    Steps: pick events with classify(); keep one event per exon with one_event_per_exon() (rbp-maps' rule: the form
    with the most inclusion reads); lift to GRCh38 with lift_event().
    A list whose content has not changed is not rewritten, so finished splice maps stay valid.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    lines = Path(se_path).read_text().splitlines()
    header = lines[0].split("\t")
    col = {name: header.index(name) for name in header}     # index() gives the first of the two ID columns
    table = pd.read_csv(se_path, sep="\t")
    assert len(table) == len(lines) - 1
    masks, lift = classify(table), lifter(cache)
    counts, written = {}, {}
    for set_name in SETS:
        picked = np.flatnonzero(masks[set_name].to_numpy())
        rows = [lines[i + 1].split("\t") for i in picked]
        keep_ids = set(one_event_per_exon(rows, col))
        kept = [r for r in rows if r[0] in keep_ids]
        lifted = [f for f in (lift_event(r, col, lift) for r in kept) if f is not None]
        out = out_dir / f"{key}.GRCh38.{set_name}.txt"
        text = lines[0] + "\n" + "".join("\t".join(f) + "\n" for f in lifted)
        if not (out.exists() and out.read_text() == text):
            out.write_text(text)
        counts[set_name] = {"rMATS events": len(rows), "one per exon": len(kept), "lifted to GRCh38": len(lifted)}
        written[set_name] = out
    return counts, written


def event_lists(key, out_dir):
    """Paths of the three GRCh38 event lists that build_event_lists() wrote for one knockdown."""
    return {s: Path(out_dir) / f"{key}.GRCh38.{s}.txt" for s in SETS}


def which_set(paths, exon):
    """The event list (set name) that holds one PICALM exon, or None."""
    c, s, a, b = EXONS[exon]["grch38"]
    for set_name, path in paths.items():
        t = pd.read_csv(path, sep="\t")
        if ((t["chr"] == c) & (t["strand"] == s) & (t["exonStart_0base"] == a) & (t["exonEnd"] == b)).any():
            return set_name
    return None


# ---------------------------------------------------------------- one exon's numbers (notebook 06)
def exon_row(se_path, exon):
    """The rMATS row of one PICALM exon (as a dict of strings).

    Matches the exon and both flanking exons. If the exon appears only with other flanks, takes the row with the
    most junction reads and says so in the 'note' field.
    """
    want = EXONS[exon]["hg19"]
    table = pd.read_csv(se_path, sep="\t", dtype=str)
    table.columns = [c if c != "ID.1" else "ID2" for c in table.columns]
    hit = table
    for k, v in want.items():
        hit = hit[hit[k] == str(v)]
    note = ""
    if len(hit) != 1:
        hit = table[(table["chr"] == want["chr"]) & (table["exonStart_0base"] == str(want["exonStart_0base"]))
                    & (table["exonEnd"] == str(want["exonEnd"]))]
        assert len(hit), f"{exon} is not in {se_path}"
        reads = hit[["IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2"]].map(lambda x: sum(_floats(x))).sum(axis=1)
        hit = hit.loc[[reads.idxmax()]]
        note = "flanking exons differ; row with the most junction reads"
    return {**hit.iloc[0].to_dict(), "note": note}


def exon_stats(row):
    """Junction counts, PSI per replicate, dPSI and FDR from an rMATS row (sample 1 = knockdown, 2 = control)."""
    ints = lambda x: [int(v) for v in _floats(x)]
    return {"IJC_kd": ints(row["IJC_SAMPLE_1"]), "SJC_kd": ints(row["SJC_SAMPLE_1"]),
            "IJC_ctl": ints(row["IJC_SAMPLE_2"]), "SJC_ctl": ints(row["SJC_SAMPLE_2"]),
            "psi_kd": _floats(row["IncLevel1"]), "psi_ctl": _floats(row["IncLevel2"]),
            "dpsi": float(row["IncLevelDifference"]), "fdr": float(row["FDR"]), "pvalue": float(row["PValue"])}


# ---------------------------------------------------------------- isoforms and sequence (notebook 07)
def gencode_transcripts(gene, cache):
    """Transcripts of one gene in GENCODE v40 (GRCh38): {transcript_id: {name, type, tags, exons}}."""
    import re
    path = fetch(GENCODE_V40, Path(cache) / "gencode.v40.annotation.gtf.gz")
    out, key = {}, 'gene_name "{}";'.format(gene)
    with gzip.open(path, "rt") as f:
        for line in f:
            if key not in line:
                continue
            p = line.rstrip("\n").split("\t")
            if p[2] not in ("transcript", "exon"):
                continue
            tid = re.search(r'transcript_id "([^"]+)"', p[8]).group(1)
            t = out.setdefault(tid, {"name": re.search(r'transcript_name "([^"]+)"', p[8]).group(1),
                                     "type": re.search(r'transcript_type "([^"]+)"', p[8]).group(1),
                                     "tags": re.findall(r'tag "([^"]+)"', p[8]), "exons": []})
            if p[2] == "exon":
                t["exons"].append([int(p[3]) - 1, int(p[4])])
    for t in out.values():
        t["exons"].sort()
    return out


def isoform_table(transcripts):
    """For each coding transcript: does it include or skip each PICALM exon?"""
    rows = []
    for tid, t in transcripts.items():
        ex = [tuple(e) for e in t["exons"]]
        row = {"transcript": tid, "name": t["name"], "type": t["type"], "MANE": "MANE_Select" in t["tags"]}
        for exon, info in EXONS.items():
            s, e = info["grch38"][2:]
            spans = ex[0][0] < s and ex[-1][1] > e
            row[info["name"] + " included"] = (s, e) in ex
            row[info["name"] + " skipped"] = spans and not any(a < e and b > s for a, b in ex)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["type", "name"]).reset_index(drop=True)


def rna_sequence(chrom, start, end, strand, cache):
    """RNA sequence of a GRCh38 interval from the UCSC API (reverse complement on the minus strand)."""
    path = fetch(f"{UCSC_API}/getData/sequence", Path(cache) / f"hg38_{chrom}_{start}_{end}.json",
                 params={"genome": "hg38", "chrom": chrom, "start": start, "end": end})
    s = json.loads(path.read_text())["dna"].upper()
    if strand == "-":
        s = s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
    return s.replace("T", "U")


# ---------------------------------------------------------------- splice maps (notebook 07)
X3 = np.arange(-300, 50)    # cassette 3' splice site region: nt from the 3'ss (negative = intron)
X5 = np.arange(-49, 301)    # cassette 5' splice site region: nt from the 5'ss (0 = last exon nt, positive = intron)


def map_png(map_dir, subject, kd_key, rep):
    """Where rbp-maps writes one splice map (its other outputs share this name)."""
    return Path(map_dir) / f"{kd_key}_KD" / f"{subject}_{rep}.se.png"


def map_outputs(map_dir, subject, kd_key, rep, lists):
    """All files a permutation-test splice map writes: the png, mean lines, normalized matrices, permutations."""
    png = map_png(map_dir, subject, kd_key, rep)
    base = str(png)[:-len(".png")]
    need = [png] + [Path(f"{base}.{lists[s].name}.{kind}") for s in SETS for kind in ("means.txt", "normed_matrix.txt")]
    need += [Path(f"{png}.{lists[SETS[2]].name}.{lists[s].name}.randsample.tsv") for s in SETS[:2]]
    return need


def map_is_finished(map_dir, subject, kd_key, rep, lists):
    """True if every output of a splice map exists and is newer than its event lists."""
    newest_list = max(p.stat().st_mtime for p in lists.values())
    return all(p.exists() and p.stat().st_mtime > newest_list for p in map_outputs(map_dir, subject, kd_key, rep, lists))


def has_read1(path, n=200_000):
    """True if a BAM still holds read 1 of paired-end eCLIP reads (rbp-maps needs read 2 only)."""
    import pysam
    with pysam.AlignmentFile(path) as bam:
        for i, read in enumerate(bam.fetch(until_eof=True)):
            if read.is_paired and read.is_read1:
                return True
            if i >= n:
                return False
    return False


def read2_only(path, threads=4):
    """Path of a read-2-only copy of a BAM (made with samtools view -f 128 if needed)."""
    if not has_read1(path):
        return str(path)
    r2 = str(path).replace(".bam", ".R2only.bam")
    if not Path(r2).exists():
        subprocess.check_call(["samtools", "view", "-@", str(threads), "-f", "128", "-b", "-o", r2, str(path)])
        subprocess.check_call(["samtools", "index", r2])
    return r2


def sense_fraction(path, events_file, n_events=3000):
    """Fraction of reads on the same strand as the cassette exon they overlap (> 0.5 means reads are sense)."""
    import pysam
    table = pd.read_csv(events_file, sep="\t").sample(frac=1, random_state=0).head(n_events)
    sense = total = 0
    with pysam.AlignmentFile(path) as bam:
        contigs = set(bam.references)
        for _, e in table.iterrows():
            if e["chr"] not in contigs:
                continue
            for read in bam.fetch(e["chr"], int(e["exonStart_0base"]), int(e["exonEnd"])):
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue
                total += 1
                sense += (read.is_reverse == (e["strand"] == "-"))
    return sense / max(total, 1), total


def _bigwigs_scaled(bam, pos, neg, sizes, direction, mapped_f260):
    """rbp-maps' own bigWig recipe (bedtools genomecov -split, RPM), with the RPM taken from the full library."""
    scale = 1e6 / float(mapped_f260)
    strands = ("-", "+") if direction == "r" else ("+", "-")       # (track for the + strand, track for the - strand)
    for strand, bw in zip(strands, (pos, neg)):
        bg = str(bw)[:-3] + ".bg"
        with open(bg, "w") as fh:
            g = subprocess.Popen(["bedtools", "genomecov", "-ibam", str(bam), "-strand", strand, "-split", "-bg",
                                  "-scale", str(scale)], stdout=subprocess.PIPE)
            s = subprocess.Popen(["sort", "-k1,1", "-k2,2n"], stdin=g.stdout, stdout=fh)
            g.stdout.close()
            assert s.wait() == 0 and g.wait() == 0, "bedtools genomecov failed for " + str(bam)
        subprocess.check_call(["bedGraphToBigWig", bg, str(sizes), str(bw)])
        os.remove(bg)


def signal_tracks(bam, background_list, signal_dir, threads=4, library=None):
    """Strand-specific RPM bigWigs for one eCLIP BAM, made with rbp-maps; returns a dict for run_splice_map().

    library: for a lite BAM (only the reads near the exons), the full BAM's read counts, {"mapped_F260": reads
    rbp-maps counts for RPM, "mapped_records": mapped records, which rbp-maps uses for its pseudocount}. With them
    the bigWigs and maps match those of the full BAM at every position the lite BAM covers.
    """
    import pysam
    from maps import plot_map
    signal_dir = Path(signal_dir)
    signal_dir.mkdir(parents=True, exist_ok=True)
    if library is None:
        bam = read2_only(bam, threads)
    frac, n = sense_fraction(bam, background_list)
    direction = "f" if frac > 0.5 else "r"
    sizes = signal_dir / f"{Path(bam).name}.chrom.sizes"
    with pysam.AlignmentFile(bam) as b:
        sizes.write_text("".join(f"{name}\t{length}\n" for name, length in zip(b.references, b.lengths)))
    pos = str(signal_dir / Path(bam).name.replace(".bam", ".norm.pos.bw"))
    neg = str(signal_dir / Path(bam).name.replace(".bam", ".norm.neg.bw"))
    out = {"bam": str(bam), "pos": pos, "neg": neg, "sense_fraction": round(frac, 3), "reads_checked": n}
    if library is None:
        plot_map.ensure_density_bigwigs(bam=bam, pos_bw=pos, neg_bw=neg, genome_file=str(sizes), direction=direction,
                                        makebigwigfiles_cmd=None, makebigwigfiles_extra_args="",
                                        makebigwigfiles_workdir=str(signal_dir / "work"))
    else:
        if not (Path(pos).exists() and Path(neg).exists()):
            _bigwigs_scaled(bam, pos, neg, sizes, direction, library["mapped_F260"])
        out["mapped"] = int(library["mapped_records"])
    return out


class _FullLibrary:
    """Stands in for the BAM handle rbp-maps keeps: it only reads .mapped (for its pseudocount)."""

    def __init__(self, mapped):
        """Stores the full library's mapped read count `mapped`."""
        self.mapped = mapped

    def close(self):
        """Does nothing: rbp-maps closes the BAM handle when a map is done, and there is no file here."""
        pass


def run_splice_map(ip, inp, lists, png, sigtest="permutation"):
    """Build one rbp-maps skipped-exon map (IP vs size-matched input) and write all its outputs next to png.

    ip, inp: dicts from signal_tracks(). Settings: IP minus input, each event scaled to sum to 1
    (per_region_subtract_and_normalize), 50 nt of exon and 300 nt of intron per region, 95% of events kept per position.
    """
    import matplotlib.pyplot as plt
    import density.ReadDensity
    import density.normalization_functions as norm
    from density import Map
    png = Path(png)
    png.parent.mkdir(parents=True, exist_ok=True)
    rd_ip = density.ReadDensity.ReadDensity(pos=ip["pos"], neg=ip["neg"], bam=ip["bam"])
    rd_in = density.ReadDensity.ReadDensity(pos=inp["pos"], neg=inp["neg"], bam=inp["bam"])
    for rd, d in ((rd_ip, ip), (rd_in, inp)):
        if "mapped" in d:                               # lite BAM: use the full library's count
            rd.bam.close()
            rd.bam = _FullLibrary(d["mapped"])
    annotations = OrderedDict((str(lists[s]), "rmats") for s in SETS)
    tests, background = [str(lists[SETS[0]]), str(lists[SETS[1]])], str(lists[SETS[2]])
    m = Map.SkippedExon(rd_ip, rd_in, str(png), norm.per_region_subtract_and_normalize, annotations,
                        exon_offset=50, intron_offset=300, min_density_threshold=0, conf=0.95)
    m.create_matrices()
    m.normalize_matrix()
    m.create_lines()
    random.seed(PERMUTATION_SEED)              # rbp-maps' background subsets: the same every run
    m.set_background_and_calculate_significance(tests, background, sigtest)
    m.write_intermediates_to_csv()
    m.plot(tests)
    plt.close("all")
    for rd in (rd_ip, rd_in):                       # each ReadDensity holds three open files; macOS allows 256
        for handle in (getattr(rd, "pos", None), getattr(rd, "neg", None), getattr(rd, "bam", None)):
            if handle is not None:
                handle.close()
    return png


def load_splice_map(map_dir, subject, kd_key, rep, lists):
    """Read a finished splice map back from disk.

    Returns dict(means={set: array}, band=(low, high), matrices={set: DataFrame}). The band is rbp-maps'
    0.5-99.5% range of 1,000 random background subsets, combined over the included and excluded sets the way
    rbp-maps combines it (lowest low, highest high).
    """
    import density.normalization_functions as norm
    png = map_png(map_dir, subject, kd_key, rep)
    base = str(png)[:-len(".png")]
    lo = hi = None
    for s in SETS[:2]:
        rs = pd.read_csv(f"{png}.{lists[SETS[2]].name}.{lists[s].name}.randsample.tsv", sep="\t", index_col=0)
        b, t = norm.bottom_top_values_from_dataframe(rs, 0.5, 0.5)
        lo = np.array(b, float) if lo is None else np.fmin(lo, np.array(b, float))
        hi = np.array(t, float) if hi is None else np.fmax(hi, np.array(t, float))
    means = {s: np.loadtxt(f"{base}.{lists[s].name}.means.txt") for s in SETS}
    matrices = {s: pd.read_csv(f"{base}.{lists[s].name}.normed_matrix.txt", index_col=0) for s in SETS}
    return {"means": means, "band": (lo, hi), "matrices": matrices}


def exon_trace(matrices, exon):
    """One PICALM exon's own row of the normalized matrix (1,400 values) and the set it is in."""
    c, s, a, b = EXONS[exon]["grch38"]
    for set_name, matrix in matrices.items():
        hits = [i for i in matrix.index if tuple(i.split("\t")[3:7]) == (c, s, str(a), str(b))]
        if hits:
            return np.asarray(matrix.loc[hits[0]], dtype=float), set_name
    return None, None


def peak(trace, region):
    """Peak position (nt) of a trace in the cassette exon's 3'ss ('3ss') or 5'ss ('5ss') region."""
    seg, x = (trace[350:700], X3) if region == "3ss" else (trace[700:1050], X5)
    return int(x[np.nanargmax(seg)])


def positions_above_band(line, band_high, region, window=(0, 100)):
    """How many of the 100 positions next to a splice site lie above the background band.

    region '3ss': the 100 intron nt before the cassette 3'ss; '5ss': the 100 intron nt after the 5'ss.
    """
    cols = slice(550, 650) if region == "3ss" else slice(750, 850)
    return int((np.asarray(line)[cols] > np.asarray(band_high)[cols]).sum())


def plot_psi(knockdown):
    """Bar chart of PSI per replicate (control vs QKI knockdown) for both PICALM exons and both cell lines."""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(EXONS), figsize=(7, 3), sharey=True)
    for ax, (exon, info) in zip(axes, EXONS.items()):
        for g, (key, kd) in enumerate(knockdown.items()):
            st = kd[exon]
            for j, (vals, color) in enumerate(((st["psi_ctl"], "#b4b3ad"), (st["psi_kd"], "#7a3db8"))):
                x = g * 2.5 + j
                ax.bar(x, np.mean(vals), width=0.8, color=color)
                ax.scatter([x - 0.12, x + 0.12][:len(vals)], vals, s=12, color="black", zorder=3)
            ax.text(g * 2.5 + 0.5, 1.02, f"ΔPSI {st['dpsi']:+.2f}", ha="center", fontsize=8)
        ax.set_xticks([g * 2.5 + 0.5 for g in range(len(knockdown))])
        ax.set_xticklabels([k.split("_")[1] + "\ncontrol | QKI KD" for k in knockdown], fontsize=8)
        ax.set_title(info["name"], fontsize=9, pad=14)
        ax.set_ylim(0, 1)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("exon inclusion (PSI)")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------- lite mode (free Colab)
def data_dirs(data):
    """Where each notebook keeps its large files. Lite mode uses its own folders, so it never mixes with a full run.

    rmats: the genome-wide rMATS results the later notebooks read (in lite mode, the full run's tables from lite_data/)
    rmats_run: where notebook 05 runs rMATS; events, maps, signal: event lists, splice maps, bigWigs.
    """
    data = Path(data)
    if LITE:
        return {"rmats": LITE_DIR / "rmats_full", "rmats_run": data / "rmats" / "lite_picalm",
                "events": data / "rbp_maps" / "event_lists_lite", "maps": data / "rbp_maps" / "maps_lite",
                "signal": data / "rbp_maps" / "bigwigs_lite"}
    return {"rmats": data / "rmats" / "rerun", "rmats_run": data / "rmats" / "rerun",
            "events": data / "rbp_maps" / "event_lists_rerun", "maps": data / "rbp_maps" / "maps_rerun",
            "signal": data / "rbp_maps" / "bigwigs"}


def lite_eclip_bams():
    """The lite eCLIP BAMs with the full libraries' read counts (lite_data/eclip/library_sizes.tsv)."""
    table = pd.read_csv(LITE_DIR / "eclip" / "library_sizes.tsv", sep="\t")
    table["bam"] = [str(LITE_DIR / "eclip" / f) for f in table["file"]]
    return table


def lite_shrna_bam(accession):
    """The lite RNA-seq BAM of one ENCODE file: only the reads on PICALM (hg19)."""
    return LITE_DIR / "shrna" / f"{accession}.PICALM.bam"
