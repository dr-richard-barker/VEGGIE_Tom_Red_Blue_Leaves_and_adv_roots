# VEGGIE Tomato (VEG-05 / OSD-767): Spaceflight Transcriptomics of *Solanum lycopersicum* cv. 'Red Robin'

Independent re-analysis of the NASA VEG-05 experiment, in which dwarf tomato
(*Solanum lycopersicum* cv. 'Red Robin') was grown aboard the International Space
Station (ISS) under **red-rich** and **blue-rich** LED light and compared with
ground controls at Kennedy Space Center. This repository holds a full
bulk RNA-seq differential-expression, functional-enrichment, cell-type
cross-reference, and light×spaceflight interaction analysis of **leaf** and
**adventitious root** tissues, together with a draft manuscript and the original
published paper for side-by-side comparison.

- **Source dataset:** NASA Open Science Data Repository [OSD-767](https://osdr.nasa.gov/bio/repo/data/studies/OSD-767) (GeneLab ID **GLDS-709**), VEG-05 experiment.
- **Reference genome:** *Solanum lycopersicum* SL3.0 / ITAG3.0.
- **Pipeline:** adapted NASA GeneLab RNAseq Consensus Pipeline (GL-DPPD-7101-G / NF_RCP 2.0.1).
- **Design:** 21 leaf samples + 15 adventitious-root samples; Flight vs Ground, with Red/Blue light as covariate.
- **Original paper:** Dixit *et al.* (2026), *BMC Plant Biology* **26**:764. DOI [10.1186/s12870-025-07621-4](https://doi.org/10.1186/s12870-025-07621-4). PDF in [`paper_original/`](paper_original/).

> **Status:** analysis complete; manuscript Results/Discussion sections in draft.
> This repository is being prepared as a citable **Zenodo** data deposit — see
> [Preparing the Zenodo deposit](#preparing-the-zenodo-deposit).

---

## Experimental design

| Tissue | Ground&nbsp;Red | Ground&nbsp;Blue | Flight&nbsp;Red | Flight&nbsp;Blue | Total |
|--------|:---:|:---:|:---:|:---:|:---:|
| Leaf   | 6 | 6 | 5 | 4 | **21** |
| Adventitious root | 3 | 3 | 5 | 4 | **15** |
| **Combined** | 9 | 9 | 10 | 8 | **36** |

Root ground controls have only 3 replicates per light condition, which limits
statistical power for detecting light×spaceflight interaction effects in root.

---

## Repository structure

```
.
├── README.md                         ← this file
├── LICENSE                           ← CC0-1.0 (public-domain dedication)
├── CITATION.cff                      ← citation metadata for Zenodo/GitHub
│
├── manuscript/                       ← draft manuscript (.docx + .md)
├── paper_original/                   ← Dixit et al. 2026 published PDF (for comparison)
├── analysis/                         ← execution-trace notebook + analysis plan
│   ├── execution_trace_analysis.ipynb
│   └── execution_trace_PLAN.md
│
└── results/
    ├── qc/                           ← MultiQC reports (raw + trimmed leaf reads)
    ├── deseq2_leaf/                  ← leaf Flight-vs-Ground DESeq2 (additive ~ light + condition)
    ├── deseq2_root/                  ← root Flight-vs-Ground DESeq2
    ├── enrichment/                   ← GO (BP/MF/CC) + KEGG enrichment, leaf & root, up & down
    ├── scrna_crossref/               ← published scRNA/atlas cell-type markers × DEGs
    ├── light_interaction/            ← light×condition interaction models (leaf/root/combined)
    ├── light_interaction_pathways/   ← KEGG pathway maps with LFC overlays (SBGN-style)
    └── cluster_annotation/           ← leaf cluster/cell-type annotation (UMAP, labels)
```

Each `results/` subfolder keeps its self-describing filename prefixes (e.g.
`root_deseq2_root_volcano.png`, `enrichment_KEGG_root_up_dotplot.png`), so files
remain unambiguous if downloaded individually from the Zenodo archive.

### File formats

| Type | Extension | Notes |
|------|-----------|-------|
| Tables | `.csv` | DEG lists, VST counts, enrichment tables, marker overlaps |
| Serialized R objects | `.rds` | DESeq2 `dds` and shrunk-result objects; load with `readRDS()` in R |
| Figures | `.png`, `.svg` | Volcano, MA, PCA, dispersion, heatmaps, dotplots, pathway maps |
| Reports | `.html`, `.md`, `.docx` | MultiQC QC, per-analysis narrative reports |
| Pathway data | `.xml` (KGML), `.json` | KEGG pathway layouts and cell-type group definitions |

A few large objects are stored zipped (`*.rds.zip`, `*.svg.zip`); unzip before use.

---

## Methods summary

1. **QC & trimming** — FastQC v0.12.1 + MultiQC v1.26; adapter/quality trimming with Trim Galore! v0.6.10 (Cutadapt v4.2), paired-end, phred33. Reports in [`results/qc/`](results/qc/).
2. **Alignment & quantification** — GeneLab RNAseq Consensus Pipeline (GL-DPPD-7101-G) against the *S. lycopersicum* SL3.0 reference; gene-level counts.
3. **Differential expression** — `DESeq2` with an **additive model `~ light + condition`** to isolate the *main effect* of spaceflight (Flight vs Ground) while adjusting for light spectrum. LFC shrinkage via `ashr`. Significance: `padj < 0.05` and `|shrunk log2FC| ≥ 1`. Outputs in [`results/deseq2_leaf/`](results/deseq2_leaf/) and [`results/deseq2_root/`](results/deseq2_root/).
4. **Functional enrichment** — `clusterProfiler` GO (BP/MF/CC, custom TERM2GENE from biomaRt) and KEGG (`organism = "sly"`); Solyc→Entrez ID mapping via biomaRt (Ensembl Plants). Outputs in [`results/enrichment/`](results/enrichment/); combined table `enrichment_enrichment_summary_all.csv`.
5. **Cell-type cross-reference** — spaceflight DEGs intersected with **published** tomato scRNA-seq / cell-atlas marker sets (leaf and root); Fisher's exact test per cell type. Outputs in [`results/scrna_crossref/`](results/scrna_crossref/). Raw 10x data was **not** reprocessed (see [Notes & limitations](#notes--limitations)).
6. **Light × spaceflight interaction** — DESeq2 models with `light:condition` term (leaf, root, and combined `~ tissue + light + condition + light:condition`); relaxed thresholds for the lower-powered interaction contrasts. Outputs in [`results/light_interaction/`](results/light_interaction/), with pathway overlays in [`results/light_interaction_pathways/`](results/light_interaction_pathways/).

---

## Key results (this re-analysis)

**Flight vs Ground, additive model (main effect of spaceflight):**

| Tissue | Significant DEGs | Up in flight | Down in flight |
|--------|:---:|:---:|:---:|
| Leaf | **2,132** | 1,032 | 1,100 |
| Adventitious root | **2,582** | 1,706 | 876 |

Adventitious roots show the larger, more up-regulation-biased response.

**Dominant enriched processes**

- **Root, up-regulated in flight:** response to oxidative stress (GO padj ≈ 6.5×10⁻¹³), hydrogen-peroxide catabolism, ethylene-activated signaling; KEGG **phenylpropanoid biosynthesis** (padj ≈ 4×10⁻³⁰), flavonoid biosynthesis, **glutathione metabolism**, stilbenoid biosynthesis.
- **Leaf, up-regulated in flight:** heme binding, monooxygenase / oxidoreductase and iron-ion binding, abscisic-acid binding; KEGG phenylpropanoid & flavonoid biosynthesis, MAPK and plant-hormone signal transduction.
- **Down-regulated:** leaf — transmembrane transport, lignin catabolism, response to water; root — RNA modification, photosynthetic electron transport (PSI).

**Cell-type localization (root, scRNA cross-reference)** — spaceflight root DEGs are most strongly enriched in **cortex**, **exodermis**, **endodermis**, **meristematic cortex**, and **xylem/vascular** marker sets (odds ratios ≈ 2.7–4.4; padj down to ≈ 7×10⁻²⁹), pointing to remodeling of the outer/ground and vascular tissues during adventitious-root formation.

**Light × spaceflight interaction** — a substantial leaf interaction signal (thousands of genes shift their flight response between Red and Blue light) but a weak root interaction (only ~15 genes at the strict threshold), consistent with the limited root ground-control replication. See `results/light_interaction/light_interaction_deseq2_summary_table.csv`.

---

## Comparison with the original paper (Dixit *et al.* 2026)

The published study and this re-analysis use the **same OSD-767 data** but **different statistical strategies**, so the raw DEG counts are not directly comparable — yet the **biology converges strongly**.

| Aspect | Dixit *et al.* 2026 (published) | This repository (re-analysis) |
|--------|--------------------------------|-------------------------------|
| DEG threshold | FDR ≤ 0.1 **and** \|logFC\| ≥ 2 | padj < 0.05 **and** \|shrunk logFC\| ≥ 1 |
| Model for spaceflight effect | Pairwise contrasts **within** each light condition, then intersected | **Additive** `~ light + condition` (light as covariate) → pooled main effect |
| Flight-vs-Ground DEGs, leaf | 45 shared across light (22 up / 23 down) | 2,132 (1,032 up / 1,100 down) |
| Flight-vs-Ground DEGs, root | 105 shared across light (24 up / 81 down) | 2,582 (1,706 up / 876 down) |
| Light-treatment shared DEGs | 198 (leaf), 305 (root) | Modeled explicitly as `light:condition` interaction |
| Roots vs leaves | Adventitious roots show "pronounced transcriptional changes" | Root has more DEGs and stronger antioxidant signal |
| Top functional themes | Oxidative-stress response, secondary-metabolite biosynthesis, hormonal regulation, cell-wall remodeling, antioxidant metabolism | Oxidative stress + H₂O₂ catabolism, phenylpropanoid/flavonoid biosynthesis, glutathione metabolism, ethylene/ABA/MAPK hormone signaling |
| Light effect | Red-rich stabilizes expression; blue-rich increases variability | Large leaf light×flight interaction; root interaction under-powered |

**Points of agreement**
- Adventitious roots are the more transcriptionally responsive tissue in both analyses.
- **Oxidative-stress / antioxidant metabolism** and **secondary-metabolite (phenylpropanoid, flavonoid) biosynthesis** are the headline enriched themes in both — the re-analysis puts precise, strongly significant GO/KEGG terms behind the paper's qualitative claims.
- **Hormone signaling** (ethylene, ABA) and **cell-wall remodeling** recur in both, and the scRNA cross-reference newly localizes the root response to cortex/exodermis/endodermis — the tissues that physically remodel during adventitious rooting.
- Light quality measurably modulates the spaceflight response in both.

**Points of difference / added value here**
- The additive-covariate model trades the paper's strict interaction-first design for **greater power on the main effect of microgravity**, which is why this analysis reports far more DEGs. Neither count is "wrong" — they answer slightly different questions (pooled main effect vs. within-light contrasts intersected).
- The looser fold-change cutoff (\|logFC\| ≥ 1 vs ≥ 2) also inflates counts relative to the paper.
- This repository adds analyses **not** in the original paper: an explicit DESeq2 **light×condition interaction** model, a published-marker **cell-type cross-reference**, and **KEGG pathway maps with expression overlays**.

> **Caveat when interpreting the comparison:** the two DEG columns above are *not* like-for-like (different thresholds, models, and contrast definitions). Treat the comparison as a **qualitative concordance check**, not a reproduction of exact gene counts.

---

## Planning, ideas & roadmap

This section captures the project plan and open ideas so collaborators (and the
Zenodo record) have the full context. It consolidates
[`analysis/execution_trace_PLAN.md`](analysis/execution_trace_PLAN.md) and the
in-analysis reports.

### Completed
- [x] QC, trimming, alignment, quantification (GeneLab pipeline, SL3.0).
- [x] DESeq2 Flight-vs-Ground for leaf and root (additive model).
- [x] GO + KEGG enrichment (up/down, both tissues).
- [x] Cell-type cross-reference of DEGs against published tomato scRNA/atlas markers.
- [x] Light×spaceflight interaction models (leaf / root / combined).
- [x] KEGG pathway visualizations with LFC overlays (SBGN-style, KGML-based).
- [x] Repository reorganization + this README.

### In progress / to do
- [ ] **Manuscript** — populate Results and Discussion sections (currently placeholders); write the Abstract Results/Conclusions; add finalized figures. See [`manuscript/`](manuscript/).
- [ ] **Root replication caveat** — foreground the 3-replicate root ground-control limitation wherever root interaction results are discussed.
- [ ] **ID-mapping coverage** — Solyc→Entrez mapping was 60–75% across gene sets; recover unmapped genes (orthology / updated ITAG4.0 mapping) to reduce enrichment blind spots.
- [ ] **Direct reproduction of the paper's contrasts** — optionally re-run with the paper's exact thresholds (FDR ≤ 0.1, \|logFC\| ≥ 2, within-light contrasts) to produce a like-for-like overlap/Venn against Dixit *et al.*
- [ ] **Provenance** — archive the analysis scripts (only the execution-trace notebook is currently included) and pin software versions in an environment file for full reproducibility.
- [ ] **Zenodo deposit** — mint DOI and cross-link to OSD-767 and the published paper.

### Ideas / stretch goals
- Reprocess raw 10x scRNA-seq for a tomato adventitious-root atlas instead of relying on published markers (impractical in the current sandbox; noted for future compute).
- Co-expression / WGCNA modules to compare against the paper's module analysis (e.g. their large ground-control-specific modules).
- Cross-species comparison of the spaceflight oxidative-stress signature with *Arabidopsis* spaceflight datasets in OSDR.

---

## Reproducibility notes

- **Load an R object:** `dds <- readRDS("results/deseq2_root/root_deseq2_root_dds.rds")`.
- **DEG tables:** `*_Flight_vs_Ground_all_genes.csv` (full results) and `*_Flight_vs_Ground_significant.csv` (thresholded) in each `deseq2_*` folder.
- **Enrichment master table:** `results/enrichment/enrichment_enrichment_summary_all.csv` (tissue × direction × database × term).
- **Interaction summary:** `results/light_interaction/light_interaction_deseq2_summary_table.csv`.

### Notes & limitations
- Raw scRNA-seq was **not** reprocessed; cell-type cross-references use published marker lists, so cell-type calls are inherited from those studies.
- `pathview`/`SBGNview` were unavailable (Rgraphviz build failure); pathway maps were generated from the KEGG REST API + KGML layouts with a custom overlay, not the standard packages.
- Analysis was run in a provisioned HPC-style sandbox (16 CPU / 64 GB RAM); some steps (marker discovery from raw 10x) were deliberately substituted with published resources.

---

## Preparing the Zenodo deposit

This repo is structured to deposit cleanly on [Zenodo](https://zenodo.org):

1. **Recommended — GitHub → Zenodo integration:** enable the repository in your
   Zenodo account, then cut a GitHub **Release** (e.g. `v1.0.0`). Zenodo archives
   the release tarball and mints a version DOI automatically. `CITATION.cff` is
   picked up for metadata.
2. **Manual upload:** create a new Zenodo upload and drag the folders in, or
   upload a zip of the working tree (exclude `.git/`). Zenodo handles the large
   `.rds`/`.csv` data files (per-file limit is generous; the whole tree is well
   within limits).
3. **Metadata to set on the record:**
   - *Upload type:* Dataset.
   - *Related identifiers:* "is derived from" OSD-767 (GLDS-709); "is supplement to" / "cites" the Dixit *et al.* 2026 DOI `10.1186/s12870-025-07621-4`.
   - *License:* CC0-1.0 (matches `LICENSE`).
   - *Keywords:* spaceflight, microgravity, tomato, *Solanum lycopersicum*, VEG-05, adventitious root, RNA-seq, DESeq2, oxidative stress, light quality.

> Before depositing, decide whether to include the `.git/` history (≈176 MB of
> binary blobs). For a clean data record, deposit the **working tree only**; keep
> the git history on GitHub.

---

## Data provenance & citation

- **Primary data:** NASA OSDR **OSD-767** / GeneLab **GLDS-709** (VEG-05). Please cite OSDR per its data-use policy.
- **Original publication:** Dixit A. R. *et al.* "Stress and light spectral quality influence the transcriptome of a tomato crop on the International Space Station." *BMC Plant Biology* (2026) **26**:764. DOI [10.1186/s12870-025-07621-4](https://doi.org/10.1186/s12870-025-07621-4).
- **This re-analysis:** see [`CITATION.cff`](CITATION.cff) (update the DOI once the Zenodo record is minted).

## License

Released under **CC0 1.0 Universal** (public-domain dedication) — see [`LICENSE`](LICENSE).
Note that the original paper PDF in `paper_original/` is © its authors under a
CC BY-NC-ND 4.0 licence and is included here only for convenience of comparison;
its reuse is governed by that licence, not CC0.
