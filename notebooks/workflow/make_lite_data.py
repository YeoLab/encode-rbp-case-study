"""Make lite_data/: the small inputs that let the workflow run on free Google Colab.

Run once, on a computer that has the full data (after a full run of notebooks 03, 04 and 05), in the
encode-rbp-case-study environment:

    python make_lite_data.py [data folder]

The data folder defaults to CASE_STUDY_DATA or <case study folder>/data. It writes
lite_data/ next to this script and prints its size (about 90 MB):

    eclip/       the seven eCLIP BAMs of the two splice maps, keeping only read 2 within 400 nt of every splice site of
                 the exons the maps use (all three groups, before one-per-exon selection) and the PICALM locus;
                 library_sizes.tsv holds each full BAM's read counts, so maps from these BAMs match the full ones
    shrna/       the eight QKI knockdown and control RNA-seq BAMs, keeping only reads on PICALM (hg19)
    rmats_full/  the full rMATS results (SE.MATS.JC.txt) of both knockdowns, from notebook 05
    gencode.v40lift37.PICALM.gtf   PICALM's gene models, so rMATS can run on the small BAMs
    regions_grch38.bed, MANIFEST.tsv (file sizes and MD5 checksums)
"""
import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import workflow_utils as U  # noqa: E402

FLANK = 400                            # rbp-maps reads 300 nt of intron and 50 nt of exon around each splice site
PICALM_GRCH38 = ("chr11", 85_950_000, 86_075_000)
PICALM_PAD = 2_000                     # around the PICALM gene in hg19, for the RNA-seq BAMs


def sh(*cmd, capture=False):
    """Run a command; return its output when capture=True.

    samtools commands fall back to the copy of samtools inside pysam when the samtools program is not installed.
    """
    cmd = [str(c) for c in cmd]
    if cmd[0] == "samtools" and not shutil.which("samtools"):
        import pysam
        args = cmd[2:]
        target = None
        if "-o" in args:                         # pysam returns the output instead of writing it
            i = args.index("-o")
            target, args = args[i + 1], args[:i] + args[i + 2:]
        out = getattr(pysam, cmd[1])(*args)
        if target:
            Path(target).write_bytes(out)
        return out.strip() if capture else None
    if capture:
        return subprocess.check_output(cmd, text=True).strip()
    subprocess.check_call(cmd)


def md5(path):
    """Returns the MD5 checksum of the file at `path`."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def splice_site_regions(se_path, cache):
    """GRCh38 intervals (BED rows) within FLANK nt of each splice site of every event in the three map groups.

    Uses all events that pass the group rules, before keeping one event per exon, so the regions cover whichever form
    the one-per-exon step keeps. Each splice site is lifted from hg19 on its own.
    """
    table = pd.read_csv(se_path, sep="\t")
    masks = U.classify(table)
    keep = masks[U.SETS[0]] | masks[U.SETS[1]] | masks[U.SETS[2]]
    lift = U.lifter(cache)
    rows = []
    for _, e in table[keep].iterrows():
        for col in ("upstreamEE", "exonStart_0base", "exonEnd", "downstreamES"):
            hit = lift.convert_coordinate(e["chr"], int(e[col]))
            if hit and len(hit) == 1:
                chrom, pos = hit[0][0], hit[0][1]
                rows.append((chrom, max(0, pos - FLANK), pos + FLANK))
    rows.append(PICALM_GRCH38)
    return sorted(set(rows))


def write_bed(rows, path):
    """Write BED rows, merging overlapping intervals."""
    merged = []
    for chrom, s, e in sorted(rows):
        if merged and merged[-1][0] == chrom and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([chrom, s, e])
    Path(path).write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in merged))
    return len(merged)


def full_eclip_bam(eclip_dir, subject, sample):
    """The full eCLIP BAM notebook 03 downloaded (not a read-2 copy made later)."""
    hits = [p for p in sorted(Path(eclip_dir).glob(f"{subject}_{sample}.genome.Aligned.sort.dedup*.bam"))
            if not p.name.endswith(".R2only.bam")]
    assert hits, f"no full BAM for {subject} {sample} in {eclip_dir}"
    return hits[0]


def lite_eclip(full, bed, out_bam):
    """Keep read 2 within the regions; return the full BAM's counts that rbp-maps uses (on its read-2 reads)."""
    paired = U.has_read1(full)
    flt = ["-f", "128"] if paired else []
    mapped_f260 = int(sh("samtools", "view", "-c", "-F", "260", *flt, full, capture=True))
    mapped_records = int(sh("samtools", "view", "-c", "-F", "4", *flt, full, capture=True))
    sh("samtools", "view", "-b", *flt, "-M", "-L", bed, "-o", out_bam, full)      # -M: jump to the regions by the index
    sh("samtools", "index", out_bam)
    kept = int(sh("samtools", "view", "-c", "-F", "260", out_bam, capture=True))
    return {"paired_input": paired, "mapped_F260": mapped_f260, "mapped_records": mapped_records, "kept_reads": kept}


def picalm_hg19(gtf):
    """PICALM's gene span in the hg19 annotation, and its GTF lines."""
    lines = [line for line in open(gtf) if 'gene_name "PICALM";' in line]
    gene = [line.split("\t") for line in lines if line.split("\t")[2] == "gene"][0]
    return (gene[0], int(gene[3]) - PICALM_PAD, int(gene[4]) + PICALM_PAD), lines


def copy_without_excluded(src, dst):
    """Copies rMATS table `src` to `dst`, leaving out the rows whose gene is in U.EXCLUDED_GENES."""
    lines = Path(src).read_text().splitlines(keepends=True)
    col = lines[0].rstrip("\n").split("\t").index("GeneID")
    keep = [lines[0]] + [l for l in lines[1:]
                         if l.split("\t")[col].strip('"').split(".")[0] not in U.EXCLUDED_GENES]
    Path(dst).write_text("".join(keep))
    print(f"rMATS {Path(src).parent.name}: {len(keep) - 1} of {len(lines) - 1} rows kept", flush=True)


def main(data):
    """Builds lite_data/ next to this script from the full data in folder `data`."""
    data = Path(data)
    out = HERE / "lite_data"
    for d in ("eclip", "shrna", "rmats_full"):
        (out / d).mkdir(parents=True, exist_ok=True)
    dirs = {"rmats": data / "rmats" / "rerun", "eclip": data / "rbp_ark_eclip_skipper2" / "full",
            "shrna": data / "downloads" / "shrna_bams_hg19", "ref": data / "downloads" / "encode_reference",
            "genome": data / "downloads" / "genome"}

    # 1. regions the splice maps read (HepG2 knockdown, the one both maps use)
    kd_keys = sorted({kd for _, kd in U.MAP_PAIRS})
    chain = [d for d in (dirs["genome"], dirs["ref"], data / "liftover")                 # a copy already on disk
             if (d / "hg19ToHg38.over.chain.gz").exists()]
    rows = []
    for kd in kd_keys:
        rows += splice_site_regions(dirs["rmats"] / kd / "SE.MATS.JC.txt", (chain or [dirs["genome"]])[0])
    bed = out / "regions_grch38.bed"
    n = write_bed(rows, bed)
    print(f"{n} regions in {bed.name}", flush=True)

    # 2. lite eCLIP BAMs
    sizes = []
    for subject in dict.fromkeys(s for s, _ in U.MAP_PAIRS):
        for sample in sorted({x for pair in U.ECLIP[subject]["reps"].values() for x in pair}):
            full = full_eclip_bam(dirs["eclip"], subject, sample)
            lite = out / "eclip" / f"{subject}_{sample}.lite.bam"
            print("eCLIP", full.name, "->", lite.name, flush=True)
            counts = lite_eclip(full, bed, lite)
            sizes.append({"subject": subject, "sample": sample, "file": lite.name, "full_bam": full.name, **counts})
    pd.DataFrame(sizes).to_csv(out / "eclip" / "library_sizes.tsv", sep="\t", index=False)

    # 3. lite RNA-seq BAMs: PICALM only
    gtf = dirs["ref"] / "gencode.v40lift37.annotation.gtf"
    (chrom, start, end), gtf_lines = picalm_hg19(gtf)
    (out / "gencode.v40lift37.PICALM.gtf").write_text("".join(gtf_lines))
    manifest_path = dirs["shrna"] / "manifest_qki.tsv"            # notebook 04; the older notebooks wrote manifest.tsv
    if not manifest_path.exists():
        manifest_path = dirs["shrna"] / "manifest.tsv"
    manifest = pd.read_csv(manifest_path, sep="\t")
    manifest = manifest[manifest["knockdown"].isin(list(U.KNOCKDOWNS))]
    assert len(manifest) == 4 * len(U.KNOCKDOWNS), manifest
    rows = []
    for _, r in manifest.iterrows():
        lite = out / "shrna" / f"{r['accession']}.PICALM.bam"
        print("RNA-seq", r["accession"], "->", lite.name, flush=True)
        bam = Path(r["bam"]) if Path(r["bam"]).exists() else dirs["shrna"] / Path(r["bam"]).name   # path from another computer
        sh("samtools", "view", "-b", "-o", lite, bam, f"{chrom}:{start}-{end}")
        sh("samtools", "index", lite)
        rows.append({**{k: r[k] for k in ("knockdown", "condition", "rep", "experiment", "accession")},
                     "file": lite.name, "reads": int(sh("samtools", "view", "-c", lite, capture=True))})
    pd.DataFrame(rows).to_csv(out / "shrna" / "manifest.tsv", sep="\t", index=False)

    # 4. the full rMATS results, without the exons of U.EXCLUDED_GENES
    for kd in U.KNOCKDOWNS:
        (out / "rmats_full" / kd).mkdir(exist_ok=True)
        copy_without_excluded(dirs["rmats"] / kd / "SE.MATS.JC.txt", out / "rmats_full" / kd / "SE.MATS.JC.txt")

    # 5. manifest
    files = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "MANIFEST.tsv")
    pd.DataFrame([{"file": str(p.relative_to(out)), "bytes": p.stat().st_size, "md5": md5(p)} for p in files]).to_csv(
        out / "MANIFEST.tsv", sep="\t", index=False)
    total = sum(p.stat().st_size for p in files)
    print(f"lite_data: {len(files)} files, {total / 1e6:.0f} MB")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CASE_STUDY_DATA", U.DATA_DEFAULT))
