# OSD-767: Leaf scRNA-seq Cluster Annotation & SBGN Pathway Projection

## Summary

This report presents two analyses for the OSD-767 spaceflight RNA-seq project:

1. **Leaf scRNA-seq cluster annotation** — 13 Leiden clusters from the GSE201931 tomato leaf cell atlas annotated with biological cell type names using marker gene characterization and cross-referencing with the published Yue et al. (2024) annotations.

2. **SBGN pathway projection** — DESeq2 shrunk log2 fold-change values for light effect (Blue vs Red) and light x flight interaction projected onto KEGG pathway maps, stratified by leaf cell types (from scRNA-seq) and root cell types (from TRAP-seq markers).

---

## 1. Leaf Cluster Annotation

### Data Source

The leaf scRNA-seq data comes from GSE201931 (Yue et al. 2024, *Plant Cell Environ* 47(7):2660-2674, DOI: 10.1111/pce.14906), which profiled 44,867 cells from healthy and tomato chlorosis virus-infected leaves using 10x Chromium. The published study identified 6 cell types: mesophyll cells, epidermal cells, guard cells, vascular cells, trichome cells, and proliferating cells.

### Annotation Method

1. **Marker gene identification**: Top 50 marker genes per cluster from Leiden clustering (CSV source, aligned with expression matrix)
2. **Gene description retrieval**: 354 unique marker genes annotated via biomaRt (Ensembl Plants), including Arabidopsis homolog names and Entrez descriptions
3. **Functional category profiling**: Computed expression of functional gene categories (photosynthesis, epidermis, defense, heat shock, etc.) across clusters
4. **KEGG pathway-level expression**: Photosynthesis pathway genes showed highest expression in Clusters 5 (mean 1.97) and 6 (mean 1.79), confirming mesophyll identity
5. **Cross-reference with Yue et al. 2024**: Mapped to the 6 published cell type categories

### Annotation Results

| Cluster | Cells | % | Sub-type | Yue et al. Type | Key Markers |
|---------|-------|---|----------|-----------------|-------------|
| 0 | 6,374 | 14.2 | Proliferating/Meristematic | Proliferating | ZF2, WRKY53, phytophthora-inhibited protease |
| 1 | 5,998 | 13.4 | Epidermis/Defense | Epidermal | EP3 endochitinase, GST, HSP90, GEM-like, dirigent |
| 2 | 4,668 | 10.4 | Epidermis (Abaxial) | Epidermal | LTP4, chitinase, peroxidase PER4, abscisic stress-ripening |
| 3 | 4,076 | 9.1 | Vascular/Phloem | Vascular | AOS, AOC2 (JA biosynthesis), hexose transporter, MAPKKK21 |
| 4 | 3,972 | 8.9 | Mesophyll (Palisade, Phenylpropanoid-active) | Mesophyll | PAL3/5, FBA1, PSBR, catalase, peroxidase |
| 5 | 3,921 | 8.7 | Mesophyll (Palisade) | Mesophyll | RBCS, PSBR, PSAO, PSAG, PSAH, RCA, PSI/PSII subunits |
| 6 | 3,566 | 8.0 | Mesophyll (Spongy) | Mesophyll | PSBR, RBCS, LHCB, PETC, PSAE, PSBP, RCA |
| 7 | 3,272 | 7.3 | Stress-responsive | Mesophyll | MBF1C (ethylene), HSP90, NQR, ClpB1, HSP17.6 |
| 8 | 3,226 | 7.2 | Stress-responsive (Heat shock) | Mesophyll | HSP70, HSP90, MBF1C, ClpB1, HSP17.6A, sucrose synthase |
| 9 | 2,219 | 4.9 | Trichome (Glandular type VI) | Trichome | acylsugar acylhydrolase, kunitz trypsin inhibitor, MYB15-like |
| 10 | 1,744 | 3.9 | Epidermis (Adaxial) | Epidermal | GDSL esterase/lipase, ACC oxidase, AOS, BURP RD22, IAA14 |
| 11 | 1,470 | 3.3 | Uncharacterized/Transitioning | Epidermal | GDSL esterase, myb-related 306, pirin, NQR, CASP-like |
| 12 | 344 | 0.8 | Uncharacterized (Small population) | Uncharacterized | metallothionein II, mostly uncharacterized genes |

### Cell Type Groupings

For pathway-level analysis, clusters were aggregated into 6 cell type groups:

| Cell Type | Clusters | Total Cells | % |
|-----------|----------|-------------|---|
| Mesophyll | 4, 5, 6, 7, 8 | 17,957 | 40.0 |
| Epidermis | 1, 2, 10, 11 | 13,880 | 30.9 |
| Proliferating | 0 | 6,374 | 14.2 |
| Vascular | 3 | 4,076 | 9.1 |
| Trichome | 9 | 2,219 | 4.9 |
| Uncharacterized | 12 | 344 | 0.8 |

**Note**: Guard cells were not detected as a separate cluster. This is expected given their small population in leaf tissue and known protoplasting bias (guard cells are resistant to enzymatic digestion and underrepresented in protoplast-based scRNA-seq).

### Key Technical Finding

A critical numbering mismatch was discovered between the JSON marker file (`leaf_cluster_markers.json`, 14 clusters 0-13) and the CSV marker file (`leaf_cluster_marker_genes.csv`, 13 clusters 0-12). The CSV markers align with the expression means matrix, while the JSON markers use a different numbering scheme. Cluster 13 in the JSON file has 0 cells in the assignments and is absent from the means matrix. All annotations were based on the CSV markers.

---

## 2. SBGN Pathway Projection

### DESeq2 Model

The bulk RNA-seq DESeq2 analysis used an interaction model: `~ light + condition + light:condition`

- **Light effect**: Blue vs Red light (main effect of light quality)
- **Condition effect**: Flight vs Ground (spaceflight effect)
- **Interaction**: Light x Condition (whether the light response differs between spaceflight and ground control)

Significance thresholds:
- Main effects: padj < 0.05, |shrunk LFC| >= 1
- Interaction: padj < 0.1, |shrunk LFC| >= 0.5

Shrunk LFC values were computed using `ashr` (adaptive shrinkage) to stabilize effect size estimates.

### DESeq2 Summary

| Model | Effect | Total Genes | Significant (strict) |
|-------|--------|-------------|---------------------|
| Leaf (interaction) | Light (B vs R) | 20,157 | 2,314 |
| Leaf (interaction) | Condition (F vs G) | 20,157 | 527 |
| Leaf (interaction) | Light x Condition | 20,157 | 2,894 (strict) / 5,042 (relaxed) |
| Root (interaction) | Light (B vs R) | 21,158 | 256 |
| Root (interaction) | Condition (F vs G) | 21,158 | 704 |
| Root (interaction) | Light x Condition | 21,158 | 15 (strict) / 44 (relaxed) |

### Pathways Selected

Six KEGG pathways were selected based on enrichment analysis results and biological relevance:

| KEGG ID | Pathway | Relevance |
|---------|---------|-----------|
| sly04075 | Plant hormone signal transduction | Core signaling; auxin, ethylene, JA pathways |
| sly00195 | Photosynthesis | Light quality directly affects photosynthetic machinery |
| sly00196 | Antenna proteins | Light-harvesting complex; strong light response expected |
| sly00940 | Phenylpropanoid biosynthesis | Secondary metabolism; stress/defense response |
| sly04712 | Circadian rhythm | Light entrainment of circadian clock |
| sly00941 | Flavonoid biosynthesis | Light-regulated secondary metabolism |

### Leaf Cell Type Pathway Diagrams

For leaf cell types, pathway gene nodes are colored by DESeq2 shrunk LFC, with node size proportional to cell-type-specific mean expression (from scRNA-seq cluster means). Four cell type groups were visualized: Mesophyll, Epidermis, Vascular, Trichome.

**12 diagrams generated** (6 pathways x 2 LFC conditions):
- `{pathway}_leaf_celltypes_light.svg/.png` — Light effect (Blue vs Red) LFC
- `{pathway}_leaf_celltypes_interaction.svg/.png` — Light x Flight interaction LFC

### Root Cell Type Pathway Diagrams

Root cell types come from the Kajala et al. (2021, *Cell* 184(12):3333-3348.e19, DOI: 10.1016/j.cell.2021.04.024) TRAP-seq translatome atlas (GSE149217), which profiled 11 cell types using INTACT/TRAP-seq. Since root cell types lack continuous expression values (only marker gene lists), pathway diagrams use:
- **Node color**: DESeq2 shrunk LFC (same as leaf)
- **Red border**: Gene is a marker for that cell type (binary indicator)
- **Node size**: Larger for marker genes, smaller for non-markers

Root cell types were grouped into 5 categories for visualization:

| Group | Constituent Cell Types | Marker Genes |
|-------|----------------------|--------------|
| Epidermis | Epidermis, Exodermis | 2,184 |
| Cortex | Cortex, Cortex_ACT2, Meristematic_Cortex | 3,115 |
| Endodermis | Endodermis | 1,396 |
| Vascular | Phloem, Vascular, Vascular_WOX4 | 2,426 |
| Xylem | Xylem | 1,167 |

Four pathways were visualized for root (excluding leaf-specific photosynthesis and antenna proteins):

**8 diagrams generated** (4 pathways x 2 LFC conditions):
- `{pathway}_root_celltypes_light.svg/.png` — Root light effect LFC
- `{pathway}_root_celltypes_interaction.svg/.png` — Root interaction LFC

### Cross-Cell-Type Comparison Heatmaps

Three summary heatmaps were generated comparing mean pathway-level LFC across all cell types:

1. **Full heatmap** (`cross_celltype_pathway_lfc_heatmap.svg`): All 6 pathways x 4 conditions (leaf light, leaf interaction, root light, root interaction) x 11 cell types (6 leaf + 5 root)

2. **Leaf-focused heatmap** (`leaf_lfc_pathway_celltype_heatmap.svg`): Leaf DESeq2 LFC conditions only, showing how leaf and root cell types respond

3. **Root-focused heatmap** (`root_lfc_pathway_celltype_heatmap.svg`): Root DESeq2 LFC conditions only

### Key Findings from Pathway Projection

**Antenna proteins (sly00196)** show the strongest light effect in leaf:
- Mean LFC ~ -2.0 across all leaf cell types (Blue light downregulates antenna genes vs Red)
- Interaction effect ~ +1.5 to +2.6 (spaceflight partially reverses the blue-light downregulation)
- Proliferating cells show the strongest interaction effect (LFC 2.63)

**Phenylpropanoid biosynthesis (sly00940)** shows strong interaction effects:
- Leaf interaction LFC ~ +1.0 to +1.3 across cell types
- Root Endodermis markers show a strong negative interaction (LFC -2.20), suggesting cell-type-specific regulation

**Circadian rhythm (sly04712)** shows cell-type-specific patterns:
- Leaf Vascular cells have the strongest light effect (LFC -1.15) and interaction (LFC +1.17)
- Root Vascular markers also show notable circadian interaction (LFC +1.76)

**Flavonoid biosynthesis (sly00941)**:
- Strong interaction effect in leaf Mesophyll (LFC +1.59) and Epidermis (LFC +1.05-1.59)
- Root cell types have minimal flavonoid pathway gene representation

**Root effects are minimal overall**:
- Root light effect LFC values are near zero across all pathways and cell types
- Root interaction effects are similarly small
- This is consistent with the DESeq2 results showing only 15-44 significant interaction genes in root

---

## 3. Technical Notes

### SBGNview Workaround

The SBGNview R package could not be installed due to an irreconcilable incompatibility between the system graphviz (v2.42, cgraph API v2) and Rgraphviz (written for cgraph API v1). As a workaround, KEGG KGML files were parsed in Python and pathway diagrams were generated using matplotlib with SBGN Process Description (PD) visual conventions:
- Gene nodes: FancyBboxPatch (rounded rectangles), colored by LFC
- Compound nodes: Circle patches
- Relations: Arrows with type-specific styling (activation = green, inhibition = red, phosphorylation = blue)

### Gene ID Mapping

DESeq2 results use ENSRNA gene IDs, while scRNA-seq data uses Solyc IDs. Mapping was performed via pre-existing entrez-annotated CSV files containing both ID types. Mapping rates: 12,627/20,157 leaf genes (62.7%) and 13,032/21,158 root genes (61.6%).

### Leaf vs Root Cell Type Data Asymmetry

Leaf cell types have continuous expression values (mean log-normalized expression per cluster from scRNA-seq), enabling quantitative node sizing. Root cell types have only binary marker gene lists (from TRAP-seq), so node sizing uses marker status (present/absent) rather than expression level. This asymmetry reflects the different technologies used in the source studies.

---

## 4. Output Files

### Cluster Annotation
- `cluster_annotation/leaf_cluster_annotations.csv` — Full annotation table (13 clusters)
- `cluster_annotation/cell_type_groups.json` — Cell type to cluster mapping
- `cluster_annotation/leaf_cluster_annotation_umap.svg/.png` — Dual-panel UMAP (clusters + cell types)

### SBGN Pathway Diagrams — Leaf
- `sbgn_pathways/sly04075_leaf_celltypes_light.svg/.png` — Hormone signaling, light LFC
- `sbgn_pathways/sly04075_leaf_celltypes_interaction.svg/.png` — Hormone signaling, interaction LFC
- `sbgn_pathways/sly00195_leaf_celltypes_light.svg/.png` — Photosynthesis, light LFC
- `sbgn_pathways/sly00195_leaf_celltypes_interaction.svg/.png` — Photosynthesis, interaction LFC
- `sbgn_pathways/sly00196_leaf_celltypes_light.svg/.png` — Antenna proteins, light LFC
- `sbgn_pathways/sly00196_leaf_celltypes_interaction.svg/.png` — Antenna proteins, interaction LFC
- `sbgn_pathways/sly00940_leaf_celltypes_light.svg/.png` — Phenylpropanoid, light LFC
- `sbgn_pathways/sly00940_leaf_celltypes_interaction.svg/.png` — Phenylpropanoid, interaction LFC
- `sbgn_pathways/sly04712_leaf_celltypes_light.svg/.png` — Circadian rhythm, light LFC
- `sbgn_pathways/sly04712_leaf_celltypes_interaction.svg/.png` — Circadian rhythm, interaction LFC
- `sbgn_pathways/sly00941_leaf_celltypes_light.svg/.png` — Flavonoid biosynthesis, light LFC
- `sbgn_pathways/sly00941_leaf_celltypes_interaction.svg/.png` — Flavonoid biosynthesis, interaction LFC

### SBGN Pathway Diagrams — Root
- `sbgn_pathways/sly04075_root_celltypes_light.svg/.png` — Hormone signaling, root light LFC
- `sbgn_pathways/sly04075_root_celltypes_interaction.svg/.png` — Hormone signaling, root interaction LFC
- `sbgn_pathways/sly00940_root_celltypes_light.svg/.png` — Phenylpropanoid, root light LFC
- `sbgn_pathways/sly00940_root_celltypes_interaction.svg/.png` — Phenylpropanoid, root interaction LFC
- `sbgn_pathways/sly04712_root_celltypes_light.svg/.png` — Circadian rhythm, root light LFC
- `sbgn_pathways/sly04712_root_celltypes_interaction.svg/.png` — Circadian rhythm, root interaction LFC
- `sbgn_pathways/sly00941_root_celltypes_light.svg/.png` — Flavonoid biosynthesis, root light LFC
- `sbgn_pathways/sly00941_root_celltypes_interaction.svg/.png` — Flavonoid biosynthesis, root interaction LFC

### Summary Heatmaps
- `sbgn_pathways/cross_celltype_pathway_lfc_heatmap.svg/.png` — Full cross-cell-type comparison
- `sbgn_pathways/leaf_lfc_pathway_celltype_heatmap.svg/.png` — Leaf LFC conditions focused
- `sbgn_pathways/root_lfc_pathway_celltype_heatmap.svg/.png` — Root LFC conditions focused

### Data Tables
- `sbgn_pathways/lfc_with_celltype_expression.csv` — Merged LFC + cell type expression (13,275 genes)

---

## 5. References

1. Yue H, Chen G, Zhang Z, et al. Single-cell transcriptome landscape elucidates the cellular and developmental responses to tomato chlorosis virus infection in tomato leaf. *Plant Cell Environ*. 2024;47(7):2660-2674. doi:10.1111/pce.14906

2. Kajala K, Gouran M, Shaar-Moshe L, et al. Innovation, conservation, and repurposing of gene function in root cell type development. *Cell*. 2021;184(12):3333-3348.e19. doi:10.1016/j.cell.2021.04.024
