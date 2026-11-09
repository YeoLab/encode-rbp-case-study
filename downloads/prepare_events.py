"""Case Study 3: knockdown event lists (GRCh38) for rbp-maps, and the regions to keep from each eCLIP BAM.

Standard library only, so it runs with any Python 3.8+.

usage: python3 prepare_events.py <rmats_dir> <hg19ToHg38.over.chain.gz> <out_dir>

<rmats_dir> holds ENCODE's QKI knockdown rMATS tars (ENCFF096LJD, ENCFF103PLX .tar.gz).
Writes <out_dir>/<RBP>_<cell>.GRCh38.<set>.txt (rMATS lines, coordinates lifted, same rules as
notebook 04) and <out_dir>/regions.bed (every splice site of every event +/- 400 nt, plus the PICALM locus).
"""
import bisect
import gzip
import io
import os
import sys
import tarfile
from pathlib import Path

KNOCKDOWNS = {"QKI_K562": "ENCFF096LJD", "QKI_HepG2": "ENCFF103PLX"}
SETS = ["included-upon-knockdown", "excluded-upon-knockdown", "native-cassette-background"]
INTERVALS = [("exonStart_0base", "exonEnd"), ("upstreamES", "upstreamEE"), ("downstreamES", "downstreamEE")]
FLANK = 400                                           # rbp-maps uses 300 nt of intron and 50 nt of exon per site
EXTRA = [("chr11", 85_960_000, 86_000_000, "PICALM")]  # locus tracks
EXCLUDED_GENES = {"ENSG00000151923"}            # exons of these genes (Ensembl IDs) are left out


# ---------- liftover (same answers as pyliftover.LiftOver.convert_coordinate) ----------
class Chain:
    def __init__(self, path):
        """Reads the gzipped UCSC chain file at `path` into alignment blocks per chromosome."""
        self.blocks = {}          # tName -> list of (tStart, tEnd, qName, qStrand, qStart, qSize, score)
        with gzip.open(path, "rt") as f:
            head = None
            for line in f:
                p = line.split()
                if not p:
                    continue
                if p[0] == "chain":
                    score, tName, tStrand, tStart = int(p[1]), p[2], p[4], int(p[5])
                    qName, qSize, qStrand, qStart = p[7], int(p[8]), p[9], int(p[10])
                    assert tStrand == "+"
                    head = [tName, tStart, qName, qSize, qStrand, qStart, score]
                    continue
                tName, t, qName, qSize, qStrand, q, score = head
                size = int(p[0])
                self.blocks.setdefault(tName, []).append((t, t + size, qName, qStrand, q, qSize, score))
                if len(p) == 3:
                    head[1], head[5] = t + size + int(p[1]), q + size + int(p[2])
        self.maxlen = {}
        for k, v in self.blocks.items():
            v.sort()
            self.maxlen[k] = max(b[1] - b[0] for b in v)
        self.starts = {k: [b[0] for b in v] for k, v in self.blocks.items()}

    def convert(self, chrom, pos):
        """Returns where 0-based position `pos` on `chrom` lifts to, as [(chrom, pos, strand, score)] like pyliftover;
        empty if it does not map.
        """
        out = []
        blocks = self.blocks.get(chrom)
        if not blocks:
            return out
        i = bisect.bisect_right(self.starts[chrom], pos) - 1
        while i >= 0 and blocks[i][0] > pos - self.maxlen[chrom] - 1:
            t0, t1, qName, qStrand, q0, qSize, score = blocks[i]
            if t0 <= pos < t1:
                q = q0 + (pos - t0)
                out.append((qName, q if qStrand == "+" else qSize - q - 1, qStrand, score))
            i -= 1
        return out


def lift_point(chain, chrom, pos):
    """Returns (chrom, position, strand) of `pos` on `chrom` lifted with `chain`, or None unless it maps to exactly one
    place.
    """
    hits = chain.convert(chrom, int(pos))
    return hits[0][:3] if len(hits) == 1 else None


def lift_fields(chain, fields, col):
    """Returns the rMATS line `fields` with its three exons lifted to GRCh38 with `chain`, or None if an exon does not
    lift whole, on the + strand of the same chromosome. `col` maps column names to indices.
    """
    chrom = fields[col["chr"]]
    out = list(fields)
    for s_name, e_name in INTERVALS:
        s, e = int(fields[col[s_name]]), int(fields[col[e_name]])
        a, b = lift_point(chain, chrom, s), lift_point(chain, chrom, e - 1)
        if a is None or b is None or a[0] != chrom or b[0] != chrom or a[2] != "+" or b[2] != "+":
            return None
        if (b[1] + 1) - a[1] != e - s:
            return None
        out[col[s_name]], out[col[e_name]] = str(a[1]), str(b[1] + 1)
    return out


# ---------- event sets ----------
def nums(text):
    """Returns the numbers in a comma-separated rMATS field `text`, skipping NA."""
    return [float(x) for x in str(text).split(",") if x not in ("NA", "")]


def classify(f, col):
    """Returns the event set of rMATS line `f` (included or excluded upon knockdown, or native background), or None;
    needs at least 30 junction reads and a gene outside EXCLUDED_GENES. `col` maps column names to indices.
    """
    if f[col["GeneID"]].strip('"').split(".")[0] in EXCLUDED_GENES:
        return None
    reads = sum(sum(nums(f[col[c]])) for c in ("IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2"))
    ctl = nums(f[col["IncLevel2"]])
    ctl_psi = sum(ctl) / len(ctl) if ctl else float("nan")
    fdr, p, d = float(f[col["FDR"]]), float(f[col["PValue"]]), float(f[col["IncLevelDifference"]])
    if reads < 30:
        return None
    if fdr < 0.05 and p < 0.05 and d >= 0.05:
        return SETS[0]
    if fdr < 0.05 and p < 0.05 and d <= -0.05:
        return SETS[1]
    if fdr >= 0.5 and 0.05 <= ctl_psi <= 0.95:
        return SETS[2]
    return None


def se_table(tar_path):
    """Returns the lines of SE.MATS.JunctionCountOnly.txt from the ENCODE rMATS tar at `tar_path`."""
    with tarfile.open(tar_path) as tar:
        for m in tar.getmembers():
            if m.name.endswith("/MATS_output/SE.MATS.JunctionCountOnly.txt"):
                return io.TextIOWrapper(tar.extractfile(m)).read().splitlines()
    raise SystemExit(f"{tar_path}: no MATS_output/SE.MATS.JunctionCountOnly.txt")


def main(rmats_dir, chain_path, out_dir):
    """Writes the GRCh38 event lists and regions.bed to `out_dir` from the rMATS tars in `rmats_dir`, lifting with the
    chain file `chain_path`.
    """
    rmats_dir, out_dir = Path(rmats_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    chain = Chain(chain_path)
    # Self-test against known coordinates (the starts of PICALM's 150-nt and 24-nt exons)
    for chrom, hg19, hg38 in () if os.environ.get("CS3_SKIP_SELFTEST") else (("chr11", 85701292, 85990249), ("chr11", 85689112, 85978069)):
        got = lift_point(chain, chrom, hg19)
        if got is None or got[1] != hg38:
            raise SystemExit(f"liftover self-test failed: {chrom}:{hg19} -> {got}, expected {hg38}")
    regions, report = [], []
    for name, acc in KNOCKDOWNS.items():
        lines = se_table(rmats_dir / f"{acc}.tar.gz")
        header = lines[0].split("\t")
        col = {}
        for i, h in enumerate(header):
            col.setdefault(h, i)            # first of the two ID columns
        by_set = {s: [] for s in SETS}
        n_hg19 = {s: 0 for s in SETS}
        for line in lines[1:]:
            f = line.split("\t")
            s = classify(f, col)
            if s is None:
                continue
            n_hg19[s] += 1
            lifted = lift_fields(chain, f, col)
            if lifted is not None:
                by_set[s].append(lifted)
        for s in SETS:
            (out_dir / f"{name}.GRCh38.{s}.txt").write_text(
                lines[0] + "\n" + "".join("\t".join(f) + "\n" for f in by_set[s]))
            report.append(f"{name}\t{s}\t{n_hg19[s]}\t{len(by_set[s])}")
            for f in by_set[s]:
                chrom = f[col["chr"]]
                for c in ("upstreamEE", "exonStart_0base", "exonEnd", "downstreamES"):
                    x = int(f[col[c]])
                    regions.append((chrom, max(0, x - FLANK), x + FLANK))
    regions += [(c, s, e) for c, s, e, _ in EXTRA]
    regions.sort()
    merged = []
    for c, s, e in regions:
        if merged and merged[-1][0] == c and s <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([c, s, e])
    (out_dir / "regions.bed").write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in merged))
    (out_dir / "event_counts.tsv").write_text("knockdown\tset\thg19\tGRCh38\n" + "\n".join(report) + "\n")
    print("knockdown\tset\thg19\tGRCh38\n" + "\n".join(report))
    print(f"regions.bed: {len(merged)} intervals, {sum(e - s for _, s, e in merged) / 1e6:.1f} Mb")


if __name__ == "__main__":
    main(*sys.argv[1:4])
