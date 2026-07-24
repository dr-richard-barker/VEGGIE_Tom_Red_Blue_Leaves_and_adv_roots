# Plan: OSD-767 — scRNA-seq Cross-Reference & Pathway Enrichment

## Background

The OSD-767 bulk RNA-seq DESeq2 analysis is **complete** for both tissues:
- **Leaf**: 2,132 significant DEGs (1,032 up, 1,100 down in Flight vs Ground)
- **Root**: 2,582 significant DEGs (1,706 up, 876 down in Flight vs Ground)

The user now requests two new analyses:
1. **Find external tomato scRNA-seq data** to identify tissue-specific (cell-type) gene clusters, then cross-reference with DEGs
2. **Run GO + KEGG pathway enrichment** on significant DEGs from each tissue

---

## Part 1: scRNA-seq Cross-Reference

### Available Datasets

After searching GEO, SRA, and literature, the most relevant published tomato scRNA-seq datasets are:

| Dataset | Tissue | Cells | Platform | Reference | GEO/Accession |
|---------|--------|-------|----------|-----------|---------------|
| Tomato leaf cell atlas | Leaf | 26,720 | 10x Chromium | Yue et al. 2024, *Plant Cell Environ* | GSE201931 |
| Tomato root cell type translatome | Root tip | N/A (INTACT/TRAP-seq, bulk per cell type) | FACS-sorted nuclei | Kajala et al. 2021, *Cell* | GSE149217 |
| Shoot-borne root scRNA-seq | Stem-borne roots | 2,671 | 10x Chromium | Omary et al. 2022, *Science* | GSE159055 |
| Root drought scRNA-seq | Seedling roots | N/A | 10x Chromium | GSE212403 | GSE212403 |
| Adventitious root scRNA-seq | Adventitious roots | 7,920 | 10x Chromium | Yang et al. 2025, *bioRxiv* | PRJCA035692/CRA022731 |

### Strategy

**Processing raw scRNA-seq data from scratch is impractical** in this sandbox (10x data requires 50+ GB per sample, alignment with STARsolo takes hours, and scanpy/Seurat analysis adds more). Instead, we will use **published cell-type marker gene lists** from these studies and cross-reference them with our DEGs.

**For Leaf (GSE201931):**
- The published paper (Yue et al. 2024) identified cell types including mesophyll, epidermis, vascular, guard cells, trichome, etc., with specific marker genes
- We will extract marker genes from the paper's supplementary data or from the processed data matrix on GEO
- Cross-reference leaf cell-type markers with leaf DEGs to identify which cell types are most affected by spaceflight

**For Root (GSE149217 — Kajala et al. 2021, Cell):**
- This is the gold-standard tomato root cell type atlas with 14+ cell types (meristem, elongation zone, cortex, endodermis, exodermis, phloem, xylem, pericycle, etc.)
- Published as a translatome atlas (INTACT/TRAP-seq), providing cell-type-enriched gene lists
- The GitHub repo (plant-plasticity/tomato-root-atlas-2020) contains processed data and marker gene lists
- Cross-reference root cell-type markers with root DEGs

**For Root (GSE159055 — Omary et al. 2022):**
- Shoot-borne root scRNA-seq with 2,671 cells across multiple cell types
- Supplementary data likely contains marker gene lists
- Secondary cross-reference if GSE149217 markers are insufficient

### Steps

1. **Fetch published marker gene lists** from supplementary tables of Yue et al. 2024 (leaf) and Kajala et al. 2021 (root)
   - Try GEO supplementary files first (processed matrices)
   - Fall back to GitHub repo for root atlas data
   - If neither available, download processed count matrix from GEO and run scanpy to identify markers ourselves

2. **Map marker gene IDs to Solyc format** (if they use different ID systems)
   - The scRNA-seq studies may use ITAG4.0/SL4.0 gene IDs vs our SL3.0/ITAG3.0
   - Use biomaRt or orthology mapping to reconcile

3. **Cross-reference markers with DEGs**
   - For each cell type, compute overlap between marker genes and DEGs
   - Fisher's exact test for enrichment of each cell type's markers in DEGs
   - Generate heatmap showing cell-type × DEG overlap significance
   - Identify which cell types are most affected by spaceflight

4. **Visualize and save**
   - Cell-type enrichment dotplot/heatmap
   - Table of DEGs annotated with cell-type identity
   - Save to `/mnt/results/scrna_crossref/`

---

## Part 2: GO + KEGG Pathway Enrichment

### Strategy

Use **clusterProfiler** (already available in base R) with biomaRt for ID mapping.

### Key Challenge: Gene ID Mapping

Our DEGs use **Solyc** gene IDs (e.g., `Solyc01g005300.3`) from Ensembl Plants SL3.0/ITAG3.0. clusterProfiler requires **Entrez Gene IDs** for `enrichGO()` and `enrichKEGG()`.

**Mapping approach:**
1. Use **biomaRt** with Ensembl Plants mart to map Solyc → Entrez Gene IDs
2. If biomaRt mapping is incomplete (common for non-model organisms), supplement with:
   - Direct KEGG API query for Solyc → Entrez mapping
   - GTF annotation parsing for genes with gene_name symbols
3. Report mapping rate; genes without Entrez IDs will be excluded from enrichment

### Steps

1. **Map Solyc gene IDs to Entrez IDs** via biomaRt
   - Connect to `plants_mart` dataset `slycopersicum_eg_gene`
   - Query: Solyc ID → Entrez Gene ID, gene_name, GO terms
   - Report mapping success rate

2. **Run GO enrichment** (BP, MF, CC) separately for each tissue × direction:
   - Leaf up-regulated DEGs (1,032 genes)
   - Leaf down-regulated DEGs (1,100 genes)
   - Root up-regulated DEGs (1,706 genes)
   - Root down-regulated DEGs (876 genes)
   - Using `clusterProfiler::enrichGO()` with:
     - `org = NULL` (custom mapping, not OrgDb)
     - `keyType = "ENTREZID"`
     - `ont = "BP"/"MF"/"CC"`
     - `pvalueCutoff = 0.05`, `qvalueCutoff = 0.05`

3. **Run KEGG pathway enrichment** separately for each tissue × direction:
   - Same 4 gene sets
   - Using `clusterProfiler::enrichKEGG()` with:
     - `organism = "sly"` (confirmed supported)
     - `keyType = "kegg"`
     - `pvalueCutoff = 0.05`, `qvalueCutoff = 0.05`

4. **Generate visualizations** (PNG format due to svglite issues):
   - GO dotplots (top 20 terms per ontology per comparison)
   - KEGG dotplots (top 20 pathways per comparison)
   - GO barplots
   - Combined GO+KEGG summary tables
   - Save to `/mnt/results/enrichment/`

5. **Export results**
   - Full enrichment tables (CSV) for each comparison
   - Summary table of top 10 terms per comparison
   - Gene-level annotation table with GO terms and KEGG pathways

---

## Execution Order

1. Set up R environment (verify clusterProfiler, biomaRt, ashr available)
2. Map Solyc → Entrez IDs via biomaRt
3. Run GO + KEGG enrichment on all 4 gene sets (leaf up, leaf down, root up, root down)
4. Generate enrichment visualizations
5. Fetch scRNA-seq marker gene lists from published studies
6. Map marker gene IDs to Solyc format
7. Cross-reference markers with DEGs, compute cell-type enrichment
8. Generate cross-reference visualizations
9. Save all results

### Compute Estimate

| Step | Time | Memory |
|------|------|--------|
| biomaRt ID mapping | 5-10 min | <1 GB |
| GO enrichment (4 sets × 3 ontologies) | 5-10 min | <1 GB |
| KEGG enrichment (4 sets) | 5-10 min | <1 GB |
| Enrichment visualizations | 5 min | <1 GB |
| scRNA-seq marker fetching | 10-30 min | <1 GB |
| Marker-DEG cross-reference | 5 min | <1 GB |
| Cross-reference visualizations | 5 min | <1 GB |
| **Total** | **~45-75 min** | **<2 GB** |

No special machine provisioning needed — default worker is sufficient for these R-based analyses.

---

## Assumptions & Caveats

1. **scRNA-seq marker genes from published studies** may use ITAG4.0 gene IDs (SL4.0 assembly), while our DEGs use ITAG3.0 (SL3.0 assembly). We will need to map between genome versions. Many Solyc IDs are stable across versions, but some differ.
2. **biomaRt mapping may be incomplete** — not all Solyc genes have Entrez IDs. We will report the mapping rate and note any bias.
3. **KEGG for tomato** uses organism code "sly" — confirmed supported by clusterProfiler. However, KEGG pathway coverage for tomato is less comprehensive than for Arabidopsis or human.
4. **GO enrichment without OrgDb**: Since there's no tomato OrgDb package, we use biomaRt for mapping and pass Entrez IDs directly to enrichGO. This works but may miss some annotations.
5. **Cell-type enrichment via Fisher's test** assumes marker gene lists are independent; some markers may be shared across cell types.
6. **GSE149217 is a translatome atlas** (INTACT/TRAP-seq), not scRNA-seq — but it provides the most comprehensive root cell-type-specific gene expression data for tomato and is the standard reference.
7. **PNG format** for plots due to svglite/systemfonts version conflict encountered previously.
