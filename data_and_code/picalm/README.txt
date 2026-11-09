Case Study 3 draft figure: QKI regulates two PICALM cassette exons in opposite directions, the 150-nt exon
(ENSE00000989276.1) and the 24-nt exon (ENSE00001416869.1).

  MAP_DIR=maps_stranded python make_fig.py      ->  CS3_PICALM_QKI_draft.png / .pdf   (numpy, pandas, matplotlib)

Panels: (a) workflow; (b) RBP-ARK windows of the 3 QKI eCLIP datasets at the 150-nt exon with GENCODE v40 isoforms;
(c) the same rows at the 24-nt exon; (d) PSI per replicate for both exons in the K562 and HepG2 QKI knockdowns, above
HepG2 junction reads for each exon; (e) rbp-maps splice maps for the HepG2 QKI knockdown with two eCLIPs (ENCORE2 K562
5007 and ENCODE3 HepG2 ENCSR570WLM, replicate 1), each with both PICALM exons' own traces.
Panel d's junction cartoons use the bar colors: control shRNA gray, QKI shRNA purple (flanking exons and the junction
with more reads).

Inputs:
  picalm_locus.json   notebook 00 (5 Nov 2026 run): RBP-ARK windows near PICALM; transcripts from GENCODE v40
  picalm.json         panel d: our rMATS 4.3.0 rerun (notebooks 02-03: Graveley-lab TopHat BAMs, GENCODE v40 GRCh37
                      mapping, fr-firststrand; junction counts, PSI per replicate, dPSI, FDR) and knockdown log2FC from
                      ENCODE DESeq2 tables without batch normalization
  picalm_encode_release.json  the same fields from ENCODE's rMATS 3.2.1 release (GENCODE v19), used by earlier drafts
  maps_stranded/      rbp-maps outputs from notebook 04 (7 Nov 2026 run): line means per set, the 0.5-99.5% background band
                      (from rbp-maps' 1,000 permutation samples), and the PICALM exon's normalized row; eCLIP replicate 1

Maps: notebook 04 on the full eCLIP BAMs, event lists from notebook 02's rMATS 4.3.0 rerun (Graveley-lab TopHat BAMs,
GENCODE v40 GRCh37 mapping, --libType fr-firststrand), filtered with FDR < 0.1, P < 0.05, |dPSI| > 0.05 (background
P > 0.5, PSI 0.05-0.95 in both conditions), one event per exon (rbp-maps subset_rmats_junctioncountonly.py), lifted to
GRCh38. Normalization per_region_subtract_and_normalize, 50 nt exon / 300 nt intron, permutation test.
