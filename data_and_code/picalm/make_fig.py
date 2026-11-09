"""Case Study 3 figure: QKI regulates two PICALM
cassette exons, the 150-nt exon (ENSE00000989276.1, Parikh et al. 2014 "exon 13") and the 24-nt exon
(ENSE00001416869.1, Parikh "exon 18").

Panels: a workflow · b 150-nt exon, eCLIP windows + isoforms · c 24-nt exon, same rows · d knockdown PSI for both exons
above the 150-nt exon's junction reads · e rbp-maps splice maps with each PICALM exon's own eCLIP trace.

Inputs (same folder unless MAP_DIR is set):
  picalm_locus.json   notebook 00: GENCODE v40 transcripts and RBP-ARK windows at PICALM (GRCh38)
  picalm.json         panel d: our rMATS 4.3.0 rerun (notebooks 02-03; junction counts, PSI per replicate, dPSI, FDR) and
                      ENCODE DESeq2 log2FC (no batch normalization); ENCODE release values in picalm_encode_release.json
  MAP_DIR             rbp-maps outputs (notebook 04 layout: <knockdown>_KD/<dataset>_<rep>.se.<knockdown>.GRCh38.<set>.txt.*)
"""
import glob
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7, "axes.titlesize": 7,
    "axes.labelsize": 7, "xtick.labelsize": 6.3, "ytick.labelsize": 6.3,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#8a8984", "axes.linewidth": 0.6,
    "xtick.color": "#52514e", "ytick.color": "#52514e",
    "text.color": "#0b0b0b", "axes.labelcolor": "#0b0b0b",
    "pdf.fonttype": 42,
})
BLUE, ORANGE, GRAY, INK, INK2, MUTED = "#2a78d6", "#eb6834", "#b4b3ad", "#0b0b0b", "#52514e", "#8a8984"
QCOL, TCOL, QDARK, TDARK = "#7a3db8", "#178a6f", "#4b2172", "#0d5242"
ODARK, BDARK = "#a8431c", "#1c5cab"
# Colors follow the rbp-maps convention: orange = more included after QKI knockdown, blue = less included.
# The 150-nt exon goes up after QKI knockdown (orange), the 24-nt exon goes down (blue).
EX150 = dict(id="ENSE00000989276.1", name="150-nt exon", key="ex13", col=ORANGE, dark=ODARK)
EX24 = dict(id="ENSE00001416869.1", name="24-nt exon", key="ex17", col=BLUE, dark=BDARK)
SEQ = LinearSegmentedColormap.from_list("seq", ["#f4f8fd", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
LOC = json.load(open("picalm_locus.json"))
KD = json.load(open("picalm.json"))["kd"]
MO = os.environ.get("MAP_DIR", "maps")
SOURCE_NOTE = os.environ.get("MAP_NOTE", "")
EX150["span"], EX24["span"] = tuple(LOC["exons"]["exon13"]), tuple(LOC["exons"]["exon24nt"])

fig = plt.figure(figsize=(7.2, 9.6))
G = gs.GridSpec(3, 2, figure=fig, height_ratios=[0.42, 0.82, 1.55], width_ratios=[1, 1.1],
                hspace=0.36, wspace=0.55, left=0.17, right=0.98, top=0.975, bottom=0.125)


def letter(x, y, s):
    """Writes panel letter `s` at figure coordinates (`x`, `y`)."""
    fig.text(x, y, s, fontsize=10, fontweight="bold", va="bottom")


# ---------- a: workflow ----------
_pa = G[0, :].get_position(fig)
axa = fig.add_axes([0.07, _pa.y0 - 0.03, 0.91, _pa.height + 0.03]); axa.axis("off"); axa.set_xlim(-1.5, 101); axa.set_ylim(-1.5, 44)
letter(0.015, _pa.y1 - 0.012, "a")


def box(x, y, w, h, title, body, edge):
    """Draws one panel-a box at (`x`, `y`), `w` wide and `h` high, with a bold `title`, the `body` text below it and an
    `edge` color.
    """
    axa.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=1.2", fc="#fcfcfb", ec=edge, lw=0.9, clip_on=False))
    axa.text(x + w / 2, y + h - 2.4, title, ha="center", va="top", fontsize=6.5, fontweight="bold")
    axa.text(x + w / 2, y + h - 7.0, body, ha="center", va="top", fontsize=5.8, color=INK2, linespacing=1.25)


def arrow(x0, y0, x1, y1):
    """Draws a panel-a arrow from (`x0`, `y0`) to (`x1`, `y1`)."""
    axa.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.8))


axa.text(0, 42, "RBP-ARK (eCLIP)", fontsize=6.6, color=BLUE, fontweight="bold")
axa.text(0, 17.2, "encodeproject.org (shRNA knockdown)", fontsize=6.6, color=ORANGE, fontweight="bold")
box(0, 21, 21, 18.5, "Gene search", "PICALM · ENCORE/ENCORE2\nQKI windows\nCoordinate View: where", BLUE)
box(25, 21, 21, 18.5, "Experiment detail", "QKI ×3 (K562 ×2, HepG2)\nQC · UCSC track hub\nreproducible windows .tsv", BLUE)
box(50, 21, 21, 18.5, "Download", "BDBag / API file list\nGlobus (bulk)", BLUE)
box(0, 0, 21, 15, "shRNA RNA-seq", "target = QKI\nK562 / HepG2 vs control", ORANGE)
box(25, 0, 21, 15, "Differential splicing", "rMATS skipped exons (SE)\nhg19 → lift to GRCh38", ORANGE)
box(50, 0, 21, 15, "ΔPSI per exon", "knockdown vs\nnon-targeting control", ORANGE)
box(77, 7, 22, 22, "Integrate", "isoforms + binding (b, c)\nsplicing change (d)\ngenome-wide rule (e)", INK)
arrow(21.6, 30, 24.4, 30); arrow(46.6, 30, 49.4, 30); arrow(21.6, 7.5, 24.4, 7.5); arrow(46.6, 7.5, 49.4, 7.5)
arrow(71.6, 30, 76.4, 23.5); arrow(71.6, 7.5, 76.4, 14.5)


# ---------- b, c: Coordinate View rows (same seven eCLIP experiments) with the isoform track, one panel per exon ----------
_pb = G[1, :].get_position(fig)
BX0, BW, CX0, CW = 0.25, 0.42, 0.715, 0.265
axb = fig.add_axes([BX0, _pb.y0, BW, _pb.height - 0.035])
axc = fig.add_axes([CX0, _pb.y0, CW, _pb.height - 0.035])
ds = [("QKI_K562_5007", "QKI · K562 · ENCORE2 5007", QDARK), ("QKI_K562_ENCSR366YOG", "QKI · K562 · ENCODE3", QDARK),
      ("QKI_HepG2_ENCSR570WLM", "QKI · HepG2 · ENCODE3", QDARK)]
WIN = {}
for subj, s, e, l2, q in LOC["windows"]:
    WIN.setdefault(subj, []).append((s, e, l2, q))
TX = LOC["transcripts"]           # GENCODE v40: {transcript_id: {name, type, tags, exons}}
ytop_w = len(ds) - 0.45
yi0 = ytop_w + 1.85
ylanes = [yi0 + 1.45, yi0]        # includes, skips


def lane(E, include):
    """Returns the protein-coding transcripts that include exon `E` (`include=True`) or skip it (`False`): their count,
    and the name, exons and MANE flag of the one drawn (MANE first, then CCDS, then by name).
    """
    picks = []
    for tid, t in TX.items():
        ex = [tuple(e) for e in t["exons"]]
        if t["type"] != "protein_coding" or not ex:
            continue
        spans = min(e[0] for e in ex) < E[0] and max(e[1] for e in ex) > E[1]
        has = E in ex
        skips = spans and not any(s < E[1] and e > E[0] for s, e in ex)
        if (include and has) or (not include and skips):
            picks.append((("MANE_Select" not in t["tags"], "CCDS" not in t["tags"], t["name"]), t["name"], ex,
                          "MANE_Select" in t["tags"]))
    picks.sort()
    n, (_, name, ex, mane) = len(picks), picks[0]
    return n, name, ex, mane


def coord_panel(ax, EX, x0, x1, key_windows, ticks, flank_text, show_labels):
    """Draws panel b or c on `ax` for exon `EX`: its including and skipping isoforms and each QKI eCLIP's windows between
    `x0` and `x1` (GRCh38). Windows in `key_windows` are labeled with their log2 enrichment; `ticks` are the x ticks,
    `flank_text` the two flank notes, and `show_labels` adds the row labels.
    """
    E = EX["span"]
    T = ax.get_yaxis_transform()
    for i, (key, lbl, col) in enumerate(ds):
        yy = len(ds) - 1 - i
        ax.add_patch(Rectangle((x0, yy - 0.32), x1 - x0, 0.64, color="#f3f2ee", lw=0, zorder=1))
        here = [w for w in WIN.get(key, []) if w[1] > x0 and w[0] < x1]
        for s, e, l2, q in here:
            shown = l2 >= 3                     # what rbp-ark.com gene search shows at its default threshold
            ax.add_patch(Rectangle((max(s, x0), yy - 0.32), min(e, x1) - max(s, x0), 0.64,
                                   color=SEQ(min(max(l2, 0), 7) / 7), lw=0, zorder=2, alpha=1 if shown else 0.55))
            if shown:
                ax.add_patch(Rectangle((s, yy - 0.32), e - s, 0.64, fill=False, ec=INK, lw=0.7, zorder=3))
            if (s, e) in key_windows and l2 >= 2:
                ax.text((s + e) / 2, yy, f"{l2:.1f}", ha="center", va="center", fontsize=4.6, zorder=4,
                        color="white" if l2 >= 3.5 else INK2)
        if not here:
            ax.text((x0 + x1) / 2, yy, "no reproducible window here", ha="center", va="center", fontsize=4.6, color=MUTED, zorder=3)
        if show_labels:
            ax.text(-0.015, yy, lbl, transform=T, ha="right", va="center", fontsize=5.4, color=col, fontweight="bold")
    for include, yl in zip((True, False), ylanes):
        n, name, exs, mane = lane(E, include)
        ax.plot([x0, x1], [yl, yl], color=MUTED, lw=0.6, zorder=2)
        step = (x1 - x0) / 5
        for xc in np.arange(x0 + step / 2, x1, step):
            if not any(s <= xc <= e for s, e in exs):
                ax.plot([xc + step * 0.06, xc, xc + step * 0.06], [yl + 0.11, yl, yl - 0.11], color=MUTED, lw=0.5, zorder=2)
        for s, e in exs:
            if e < x0 or s > x1:
                continue
            ax.add_patch(Rectangle((s, yl - 0.31), e - s, 0.62, color=EX["col"] if (s, e) == E else GRAY, lw=0, zorder=3))
        ax.text(x1 - (x1 - x0) * 0.015, yl + 0.5, f"{n} coding ({'MANE ' if mane else 'e.g. '}{name})", fontsize=4.6,
                color=INK2, ha="left", va="center")
    ax.text(x1 - (x1 - x0) * 0.015, ylanes[1] - 0.72, flank_text[0], fontsize=4.6, color=MUTED, ha="left", va="center")
    ax.text(x0 + (x1 - x0) * 0.015, ylanes[1] - 0.72, flank_text[1], fontsize=4.6, color=MUTED, ha="right", va="center")
    ytop = ylanes[0] + 0.62
    ax.add_patch(Rectangle((E[0], -0.5), E[1] - E[0], ytop + 0.5, color=EX["col"], alpha=0.12, lw=0, zorder=0))
    ax.text((E[0] + E[1]) / 2, ytop + 0.3, f"{EX['name']}\n({EX['id']})", ha="center", va="bottom", fontsize=5.4,
            color=EX["dark"], fontweight="bold", linespacing=1.15)
    ax.set_xlim(x1, x0); ax.set_ylim(-0.7, ytop + 1.9)
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xticks(ticks); ax.set_xticklabels([f"{t:,}" for t in ticks])
    ax.tick_params(axis="x", labelsize=5.6)


# 150-nt exon: 0-100 nt windows upstream of the 3'ss and downstream of the 5'ss; 24-nt exon: 0-100 nt downstream of the 5'ss
coord_panel(axb, EX150, 85989800, 85990850, {(85990399, 85990499), (85990149, 85990249)},
            [85990000, 85990250, 85990500, 85990750], ("← exon 12: 6.4 kb", "exon 14: 6.3 kb →"), True)
coord_panel(axc, EX24, 85977750, 85978400, {(85977969, 85978068)},
            [85977800, 85978000, 85978200, 85978400], ("← upstream exon: 3.0 kb", "1.4 kb →"), False)
TB = axb.get_yaxis_transform()
axb.text(-0.015, ytop_w + 0.1, "eCLIP windows", transform=TB, ha="right", va="center", fontsize=5.4, color=INK2, fontweight="bold")
axb.text(-0.015, ylanes[0] + 0.12, "includes the exon", transform=TB, ha="right", va="center", fontsize=5.4, color=INK, fontweight="bold")
axb.text(-0.015, ylanes[1] + 0.12, "skips the exon", transform=TB, ha="right", va="center", fontsize=5.4, color=INK, fontweight="bold")
axb.text(-0.015, ylanes[0] + 1.05, "GENCODE v40 isoforms (− strand)", transform=TB, ha="right", va="center", fontsize=5.4,
         color=INK2, fontweight="bold")
n_short = sum(1 for t in TX.values() if t["type"] == "protein_coding" and [EX150["span"][0], EX150["span"][1] - 21] in t["exons"])
if n_short:
    axb.text(85989800 + 15, ylanes[0] + 0.5, f"+{n_short} use a 129-nt form (alt. 3′ss)", fontsize=4.6, color=MUTED,
             ha="right", va="center")
axb.set_xlabel("chr11 (GRCh38), drawn 5′ → 3′ on the minus strand", fontsize=6.2)
axc.set_xlabel("chr11 (GRCh38), 5′ → 3′, 11.6 kb downstream of b", fontsize=6.2)
_tb = axb.get_position(); _tc = axc.get_position()
letter(0.015, _tb.y1 + 0.012, "b"); letter(_tc.x0 - 0.03, _tc.y1 + 0.012, "c")
fig.text(_tb.x0 - 0.2, _tb.y1 + 0.022, "150-nt exon: QKI binds upstream", fontsize=6.8, fontweight="bold",
         color=ODARK, va="bottom")
fig.text(_tb.x0 - 0.2, _tb.y1 + 0.008, "RBP-ARK reproducible windows (Skipper 2.0.0), all three QKI eCLIPs · "
         "outlined = shown by gene search (log2 enrichment ≥ 3)", fontsize=5.0, color=INK2, va="bottom")
fig.text(_tc.x0, _tb.y1 + 0.022, "24-nt exon: QKI binds downstream", fontsize=6.8, fontweight="bold", color=BDARK, va="bottom")
cax = fig.add_axes([0.085, _tb.y0 - 0.032, 0.13, 0.006])
cb = fig.colorbar(plt.cm.ScalarMappable(cmap=SEQ, norm=plt.Normalize(0, 7)), cax=cax, orientation="horizontal")
cb.set_label("window log2 enrichment", fontsize=4.8, labelpad=1); cb.ax.xaxis.set_label_position("top")
cb.ax.tick_params(labelsize=4.6, length=2, pad=1); cb.outline.set_visible(False)


# ---------- d: knockdown PSI for both exons (top) and each exon's junction reads in HepG2 (bottom) ----------
_pd = G[2, 0].get_position(fig)
axd = fig.add_axes([0.075, _pd.y0 + _pd.height * 0.56, _pd.x1 - 0.075, _pd.height * 0.33])
groups = [("QKI|K562", "K562"), ("QKI|HepG2", "HepG2")]
BLOCK = 5.6
for b, EX in enumerate((EX150, EX24)):
    xb = b * BLOCK
    for g, (key, cell) in enumerate(groups):
        e = KD[key][EX["key"]]
        xg = xb + g * 2.3
        for j, (vals, c) in enumerate(((e["psi_ctl"], GRAY), (e["psi_kd"], QCOL))):
            axd.bar(xg + j * 0.9, np.mean(vals), width=0.8, color=c, zorder=2, lw=0)
            axd.scatter([xg + j * 0.9 - 0.14, xg + j * 0.9 + 0.14], vals, s=5, color=INK, zorder=3, lw=0)
        fdr = "n.s." if e["fdr"] >= 0.05 else ("FDR <1e-10" if e["fdr"] < 1e-10 else f"FDR {e['fdr']:.0e}".replace("e-0", "e-"))
        axd.text(xg + 0.45, 1.0, f"ΔPSI {e['dpsi']:+.2f}\n{fdr}", ha="center", va="bottom", fontsize=4.6, color=INK2,
                 linespacing=1.1)
        axd.text(xg, -0.04, "ctl", ha="center", va="top", fontsize=4.5, color=INK2, transform=axd.get_xaxis_transform())
        axd.text(xg + 0.9, -0.04, "KD", ha="center", va="top", fontsize=4.5, color=QCOL, transform=axd.get_xaxis_transform())
        axd.text(xg + 0.45, -0.15, f"{cell}\nlog2FC {KD[key]['log2fc']:.2f}", ha="center", va="top", fontsize=4.7,
                 color=INK if cell == "HepG2" else INK2, fontweight="bold" if cell == "HepG2" else None,
                 transform=axd.get_xaxis_transform(), linespacing=1.15)
    axd.text(xb + 1.6, 1.27, f"{EX['name']} ({'b' if EX is EX150 else 'c'})", ha="center", va="bottom", fontsize=5.8,
             color=EX["dark"], fontweight="bold", transform=axd.get_xaxis_transform())
    axd.plot([xb - 0.5, xb + 3.7], [1.255, 1.255], color=EX["col"], lw=1.4, transform=axd.get_xaxis_transform(), clip_on=False)
axd.set_xticks([]); axd.set_ylim(0, 1); axd.set_xlim(-0.65, BLOCK + 3.75); axd.set_yticks([0, 0.5, 1])
axd.set_ylabel("exon inclusion (PSI)", fontsize=6.2)
axd.axvline(BLOCK - 0.95, color="#e3e2de", lw=0.6)
fig.text(0.075 + (_pd.x1 - 0.075) / 2, _pd.y0 + _pd.height * 0.56 - 0.05,
         "QKI shRNA vs control shRNA; rMATS ΔPSI and FDR, 2 replicates each (dots); log2FC = QKI mRNA (DESeq2)",
         ha="center", va="top", fontsize=4.6, color=MUTED)
letter(0.015, _pd.y1 + 0.022, "d")
fig.text(0.075, _pd.y1 + 0.022, "QKI knockdown changes both exons, in opposite directions", fontsize=6.8, fontweight="bold",
         va="bottom")

axs = fig.add_axes([0.075, _pd.y0 - 0.005, _pd.x1 - 0.075, _pd.height * 0.36]); axs.axis("off")
U, C, Dn = (0, 18), (44, 54), (80, 98)
CART = [(EX150, ("exon 12", "exon 14"), 5.15), (EX24, ("exon 17", "exon 19"), 1.15)]
tot = max(max(sum(KD["QKI|HepG2"][EX["key"]][k]) for k in ("IJC_ctl", "SJC_ctl", "IJC_kd", "SJC_kd")) for EX, _, _ in CART)


def arc(xa, xb, y, h, lw, col, up):
    """Draws one panel-d junction arc from `xa` to `xb` at height `y`, rising `h` above the line (`up`) or below it, with
    line width `lw` and color `col`.
    """
    xs = np.linspace(xa, xb, 40); ys = y + (1 if up else -1) * h * np.sin(np.linspace(0, np.pi, 40))
    axs.plot(xs, ys, color=col, lw=lw, solid_capstyle="round")


for EX, (lu, ld), ytop_c in CART:
    e = KD["QKI|HepG2"][EX["key"]]
    for name, inc, skp, y, col in [("control shRNA", sum(e["IJC_ctl"]), sum(e["SJC_ctl"]), ytop_c - 1.05, INK),
                                   ("QKI shRNA", sum(e["IJC_kd"]), sum(e["SJC_kd"]), ytop_c - 2.25, QCOL)]:
        axs.plot([U[0], Dn[1]], [y, y], color=INK2, lw=0.6)
        # colors match the bars above: control = gray, QKI shRNA = purple (flanking exons and the junction with more
        # reads); the junction with fewer reads is a thin dark line; in the QKI row the cassette keeps its exon color
        cond = QCOL if name == "QKI shRNA" else GRAY
        for (s0, e0), c in ((U, cond), (C, EX["col"] if name == "QKI shRNA" else GRAY), (Dn, cond)):
            axs.add_patch(Rectangle((s0, y - 0.1), e0 - s0, 0.2, color=c, lw=0, zorder=3))
        lwi, lws = 0.4 + 3.2 * inc / tot, 0.4 + 3.2 * skp / tot
        inc_c, skp_c = (cond, INK2) if inc >= skp else (INK2, cond)
        arc(U[1], C[0], y + 0.1, 0.32, lwi, inc_c, True); arc(C[1], Dn[0], y + 0.1, 0.32, lwi, inc_c, True)
        arc(U[1], Dn[0], y - 0.1, 0.36, lws, skp_c, False)
        axs.text(Dn[1] + 4, y + 0.17, f"incl. {inc}", va="center", fontsize=5.2, color=EX["dark"], fontweight="bold")
        axs.text(Dn[1] + 4, y - 0.2, f"skip {skp}", va="center", fontsize=5.2, color=INK2)
        axs.text(U[0] - 3, y, name, ha="right", va="center", fontsize=5.4, color=col,
                 fontweight="bold" if name == "QKI shRNA" else None)
    axs.text(sum(U) / 2, ytop_c - 0.3, lu, ha="center", fontsize=5.2, color=INK2)
    axs.text(sum(C) / 2, ytop_c - 0.3, EX["name"], ha="center", fontsize=5.2, color=EX["dark"], fontweight="bold")
    axs.text(sum(Dn) / 2, ytop_c - 0.3, ld, ha="center", fontsize=5.2, color=INK2)
axs.plot([U[0] - 30, Dn[1] + 22], [2.2, 2.2], color="#e3e2de", lw=0.6, clip_on=False)
axs.text(Dn[1] + 4, CART[0][2] - 0.3, "reads", fontsize=5.0, color=MUTED)
axs.set_xlim(-30, 116); axs.set_ylim(-1.7, 5.9)
axs.text(-30, 5.9, "Junction reads in HepG2 (2 replicates summed; not to scale)", fontsize=5.6, color=INK2, va="top")


# ---------- e: rbp-maps splice maps (QKI) with each PICALM exon's own trace ----------
MAPS = [("QKI", "QKI_K562_5007", "QKI_HepG2", "rep1", "QKI ENCORE2\n(K562) eCLIP ×\nHepG2 KD", QDARK,
         [(1, "exon13", EX150), (2, "exon24nt", EX24)]),
        ("QKI", "QKI_HepG2_ENCSR570WLM", "QKI_HepG2", "rep1", "QKI ENCODE3\n(HepG2) eCLIP ×\nHepG2 KD", QDARK,
         [(1, "exon13", EX150), (2, "exon24nt", EX24)])]
SETS = ["included-upon-knockdown", "excluded-upon-knockdown", "native-cassette-background"]
TRACE_COL = {("QKI", "exon13"): (ORANGE, ODARK), ("QKI", "exon24nt"): (BLUE, BDARK)}


def kd_dir(kd):
    """Returns the splice-map folder of knockdown `kd` under MAP_DIR (`<kd>_KD`)."""
    return next((d for d in (f"{MO}/{kd}_KD", f"{MO}/{kd}KD") if os.path.isdir(d)), f"{MO}/{kd}_KD")


def path(ds_, kd, rep, s, kind):
    """Returns the path of one rbp-maps output: eCLIP dataset `ds_`, knockdown `kd`, replicate `rep`, event set `s` and
    file `kind` (for example `means.txt`).
    """
    return f"{kd_dir(kd)}/{ds_}_{rep}.se.{kd}.GRCh38.{s}.txt.{kind}"


def picalm_trace(ds_, kd, rep, exon):
    """Returns the PICALM exon's own normalized eCLIP signal (1,400 values) for dataset `ds_`, knockdown `kd`, replicate
    `rep` and `exon`: notebook 04's saved trace, or else its row of the normalized matrix (None if absent).
    """
    saved = f"{kd_dir(kd)}/PICALM_{exon}_{ds_}_{rep}.txt"       # written by notebook 04
    if os.path.exists(saved):
        return np.loadtxt(saved)
    s0, e0 = LOC["exons"][exon]
    for s in SETS:
        m = pd.read_csv(path(ds_, kd, rep, s, "normed_matrix.txt"), index_col=0)
        hit = [i for i in m.index if f"\t{s0}\t{e0}\t" in i]
        if hit:
            return m.loc[hit[0]].to_numpy(float)
    return None


def band(ds_, kd, rep, n_perm=1000, seed=0):
    """0.5-99.5% band of background means over random subsets the size of each test set (rbp-maps' permutation
    band). Uses rbp-maps' randsample files when the permutation test was run, else samples the background matrix."""
    out = f"{kd_dir(kd)}/band_{ds_}_{rep}.tsv"
    if os.path.exists(out):
        d = np.loadtxt(out); return d[:, 0], d[:, 1]
    lo = np.full(1400, np.inf); hi = np.full(1400, -np.inf)
    rs = sorted(glob.glob(f"{kd_dir(kd)}/{ds_}_{rep}.se.png.*randsample.tsv"))
    if rs:
        import density.normalization_functions as rnorm      # from rbp-maps (github.com/YeoLab/rbp-maps)
        for f in rs:
            b, t = rnorm.bottom_top_values_from_dataframe(pd.read_csv(f, sep="\t", index_col=0), 0.5, 0.5)
            lo = np.minimum(lo, np.array(b, float)); hi = np.maximum(hi, np.array(t, float))
    else:
        bg = pd.read_csv(path(ds_, kd, rep, SETS[2], "normed_matrix.txt"), index_col=0).to_numpy(float)
        rng = np.random.default_rng(seed)
        for s in SETS[:2]:
            n = len(pd.read_csv(path(ds_, kd, rep, s, "normed_matrix.txt"), index_col=0, usecols=[0]))
            means = np.array([np.nanmean(bg[rng.choice(len(bg), n, replace=False)], axis=0) for _ in range(n_perm)])
            lo = np.minimum(lo, np.nanpercentile(means, 0.5, axis=0)); hi = np.maximum(hi, np.nanpercentile(means, 99.5, axis=0))
    np.savetxt(out, np.column_stack([lo, hi]))
    return lo, hi


sub = G[2, 1].subgridspec(4, 4, hspace=0.32, wspace=0.08, height_ratios=[1, 0.46, 1, 0.46])
# positions: intron nucleotides are numbered from 1 away from the splice site (-1 = last intron nt before a 3′ss,
# +1 = first intron nt after a 5′ss), as in notebook 04
regions = [("upstream\nexon end", -49, "ex_left"), ("cassette\n3′ss", -300, "ex_right"), ("cassette\n5′ss", -49, "ex_left"),
           ("downstream\nexon 3′ss", -300, "ex_right")]
SC = 1e4
axes_e = []
for k, (rbp, ds_, kd, rep, ylab, dark, traces) in enumerate(MAPS):
    inc = np.loadtxt(path(ds_, kd, rep, SETS[0], "means.txt")) * SC
    exc = np.loadtxt(path(ds_, kd, rep, SETS[1], "means.txt")) * SC
    bgm = np.loadtxt(path(ds_, kd, rep, SETS[2], "means.txt")) * SC
    lo, hi = [v * SC for v in band(ds_, kd, rep)]
    ymax = np.nanmax([inc.max(), exc.max(), np.nanmax(hi)]) * 1.25
    ymin = min(0, np.nanmin([inc.min(), exc.min(), np.nanmin(lo)]) * 1.1)
    tr = {}
    for j, exon, EX in traces:
        t_ = picalm_trace(ds_, kd, rep, exon) * SC
        seg = slice(350, 700) if j == 1 else slice(700, 1050)
        tx = np.arange(350) + (-300 if j == 1 else -49)
        tr[j] = (t_[seg], tx, int(tx[np.nanargmax(t_[seg])]), exon, EX)
    for j, (lab, x0r, exside) in enumerate(regions):
        ax = fig.add_subplot(sub[2 * k, j]); axes_e.append(ax)
        sl = slice(j * 350, (j + 1) * 350); x = np.arange(350) + x0r
        if exside == "ex_left": ax.axvspan(-49.5, 0.5, color="#f3f2ee", lw=0, zorder=0)
        else: ax.axvspan(0, 50, color="#f3f2ee", lw=0, zorder=0)
        ax.fill_between(x, lo[sl], hi[sl], color=GRAY, alpha=0.45, lw=0, zorder=1)
        ax.plot(x, bgm[sl], color=INK, lw=0.6, zorder=2)
        ax.plot(x, exc[sl], color=BLUE, lw=0.8, zorder=3)
        ax.plot(x, inc[sl], color=ORANGE, lw=0.9, zorder=4)
        for arr, c, yy in [(inc, ORANGE, 0.97), (exc, BLUE, 0.91)]:
            sig = arr[sl] > hi[sl]
            ax.scatter(x[sig], np.full(sig.sum(), ymin + (ymax - ymin) * yy), s=1.2, color=c, marker="s", lw=0)
        ax.set_ylim(ymin, ymax); ax.set_xlim(x[0], x[-1])
        ax.set_xticks([0]); ax.set_xticklabels([])
        if j in tr:
            ticks_ = [-300, -200, -100, 0] if j == 1 else [0, 100, 200, 300]
            ax.set_xticks(ticks_); ax.set_xticklabels([])
            for xv in ticks_:
                ax.axvline(xv, color="#e3e2de", lw=0.5, zorder=0.5)
            ax.axvline(tr[j][2], color=TRACE_COL[(rbp, tr[j][3])][1], lw=0.7, ls=(0, (2, 2)), zorder=5)
        ax.tick_params(labelsize=5.2, length=2)
        if k == 0:
            ax.set_title(lab, fontsize=5.6, pad=2)
        if j > 0:
            ax.set_yticks([]); ax.spines["left"].set_visible(False)
        else:
            ax.set_ylabel(f"{ylab}\nIP−input (×10⁻⁴)", fontsize=4.9, color=dark)
            ax.tick_params(axis="y", labelsize=5.2)
    for j in range(4):
        axt = fig.add_subplot(sub[2 * k + 1, j])
        if j not in tr:
            axt.axis("off")
            # label the trace from the empty neighbouring cell: left of a 3'ss trace, right of a 5'ss trace
            for jt, (_, _, _, exon, EX) in tr.items():
                tdark = TRACE_COL[(rbp, exon)][1]
                if (jt == 1 and j == 0) or (jt == 2 and j == 1 and 1 not in tr):
                    axt.text(0.98, 0.5, f"PICALM {EX['name']}\n({'b' if EX is EX150 else 'c'}, d) →", transform=axt.transAxes,
                             ha="right", va="center", fontsize=4.9, color=tdark, fontweight="bold")
                if jt == 2 and j == 3 and 1 in tr:
                    axt.text(0.04, 0.5, f"← PICALM\n{EX['name']}\n({'b' if EX is EX150 else 'c'}, d)", transform=axt.transAxes,
                             ha="left", va="center", fontsize=4.9, color=tdark, fontweight="bold")
            continue
        tt, tx, pk, exon, EX = tr[j]
        tcol, tdark = TRACE_COL[(rbp, exon)]
        if j == 1: axt.axvspan(0, 50, color="#f3f2ee", lw=0, zorder=0)
        else: axt.axvspan(-49.5, 0.5, color="#f3f2ee", lw=0, zorder=0)
        axt.fill_between(tx, 0, np.clip(tt, 0, None), color=tcol, alpha=0.35, lw=0)
        axt.plot(tx, tt, color=tdark, lw=0.8)
        axt.axhline(0, color=MUTED, lw=0.4)
        axt.set_xlim(tx[0], tx[-1])
        ticks_ = [-300, -200, -100, 0] if j == 1 else [0, 100, 200, 300]
        axt.set_xticks(ticks_)
        axt.set_xticklabels(["", "−200", "", "0"] if j == 1 else ["0", "", "200", ""])
        for xv in ticks_:
            axt.axvline(xv, color="#e3e2de", lw=0.5, zorder=0.5)
        axt.axvline(pk, color=tdark, lw=0.7, ls=(0, (2, 2)), zorder=5)
        axt.text(0.03 if j == 1 else 0.97, 0.92, f"peak {'−' if pk < 0 else '+'}{abs(pk)} nt\n(max {np.nanmax(tt):.0f})",
                 transform=axt.transAxes, fontsize=4.5, color=tdark, va="top", ha="left" if j == 1 else "right", linespacing=1.1)
        axt.tick_params(labelsize=5.0, length=2)
        axt.set_yticks([])
_pe = G[2, 1].get_position(fig)
letter(_pe.x0 - 0.115, _pd.y1 + 0.022, "e")
fig.text(_pe.x0 - 0.085, _pd.y1 + 0.022, "rbp-maps splice maps and the PICALM exons' own traces", fontsize=6.8, fontweight="bold",
         va="bottom")
fig.legend(handles=[Line2D([], [], color=ORANGE, lw=1.2), Line2D([], [], color=BLUE, lw=1.2),
                          Rectangle((0, 0), 1, 1, color=GRAY, alpha=0.45)],
                 labels=["more included after KD", "less included after KD", "background 0.5–99.5%"],
                 loc="upper center", bbox_to_anchor=(_pe.x0 + _pe.width / 2 - 0.03, 0.095), ncol=2, fontsize=5.0, frameon=False, handlelength=1.2,
                 columnspacing=0.8)
fig.text(_pe.x0 + _pe.width / 2 - 0.03, 0.068, "50 nt exon + 300 nt intron around each splice site.\n"
         "Each PICALM trace (signal at that one exon, same units) is colored\n"
         "by its exon's change after that knockdown. eCLIP replicate 1 shown;\n"
         "replicate 2 maps are saved alongside replicate 1.",
         fontsize=4.8, ha="center", va="top", linespacing=1.3)
fig.text(0.5, 0.01, "DRAFT Case Study 3 figure · data from rbp-ark.com (Skipper 2.0.0) and encodeproject.org"
         + ("\n" + SOURCE_NOTE if SOURCE_NOTE else ""), ha="center", fontsize=5.6, color=MUTED, linespacing=1.4)
out = "CS3_PICALM_QKI_draft"
fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=220)
print("saved", out)
