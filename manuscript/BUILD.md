# Manuscript build notes

Three renderings of the same manuscript live here:

| File | Format | Notes |
|------|--------|-------|
| `Transcriptomic Profiling ... Covariates.md` | Markdown | Canonical source; figures embedded as base64 |
| `Transcriptomic Profiling ... Covariates.docx` | Word | Regenerated from the `.md` with pandoc |
| `manuscript_npj.tex` | LaTeX | npj Microgravity–style, uses high-resolution figures in `figures/` |

## Regenerate the DOCX from the Markdown
```bash
pandoc "Transcriptomic Profiling of Dwarf Tomato (Solanum lycopersicum cv. Red Robin) in Spaceflight_ Evaluating the Main Effects of Microgravity Under Red and Blue Light Covariates.md" \
  -o  "Transcriptomic Profiling of Dwarf Tomato (Solanum lycopersicum cv. Red Robin) in Spaceflight_ Evaluating the Main Effects of Microgravity Under Red and Blue Light Covariates.docx"
```
(The base64 images in the `.md` are embedded automatically.)

## Build the LaTeX PDF
```bash
pdflatex manuscript_npj.tex
pdflatex manuscript_npj.tex   # second pass resolves \ref and citations
```
Requires a standard TeX distribution (TeX Live / MacTeX). The preamble uses only
common packages (`geometry, helvet, graphicx, subcaption, caption, booktabs,
authblk, titlesec, natbib, hyperref`). Figures are read from `figures/` — raster
panels as PNG, leaf DESeq2 panels as PDF (vector, converted from the repo SVGs).

## Figure sources (`figures/`)
The LaTeX figures are the high-resolution originals from `../results/`, not the
low-resolution thumbnails embedded in the `.md`/`.docx`:

- Fig 1 — leaf/root PCA + sample-distance heatmaps (`deseq2_leaf`, `deseq2_root`)
- Fig 2 — leaf volcano + MA (`deseq2_leaf`, SVG→PDF)
- Fig 3 — root volcano + MA (`deseq2_root`)
- Fig 4 — root GO/KEGG enrichment dot plots (`enrichment`)
- Fig 5 — leaf GO/KEGG enrichment dot plots (`enrichment`)
- Fig 6 — root DEG × cell-type cross-reference (`scrna_crossref`)

To submit to npj Microgravity, the body text and figure/table blocks can be moved
into the official Springer Nature class (`sn-jnl.cls`) with minimal edits; the
current preamble emulates that house style so the draft is readable as-is.
