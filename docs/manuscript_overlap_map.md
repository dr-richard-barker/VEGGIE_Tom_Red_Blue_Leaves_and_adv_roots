# Manuscript overlap map — VEG-05 / OSD-767

**Internal planning document — not for submission.** This maps the overlap between the
two manuscripts derived from the same OSD-767 dataset, to support the one-paper-vs-two
decision. It deliberately lives outside either submission package so a reviewer never
sees a self-overlap assessment.

- **Manuscript A (this repo):** *Transcriptomic profiling of dwarf tomato in spaceflight:
  main effects of microgravity under red- and blue-light covariates* — `manuscript/`.
- **Manuscript B (sibling repo):** *Light quality modulates spaceflight stress decoding
  with cell-type asymmetry in tomato* —
  [`PhysioSpace_stress_decoding_VEG05`](https://github.com/dr-richard-barker/PhysioSpace_stress_decoding_VEG05)
  (formerly `DeepLearning_VEG05`).

Prepared 2026-07-24.

---

## Shared foundation (identical lineage)

Both manuscripts are built on exactly the same upstream analysis:

- Same primary data: NASA OSD-767 / GLDS-709 (VEG-05), leaf (n=21) + adventitious root (n=15).
- Same pipeline: GeneLab RNAseq Consensus Pipeline → RSEM → DESeq2 → clusterProfiler (GO/KEGG) → Fisher's exact cell-type cross-reference.
- Same cell-type markers: Yue et al. (leaf scRNA) + Kajala et al. (root TRAP-seq).
- **Byte-identical shared files:** `execution_trace_PLAN.md` and the entire `scDATA/` /
  `scrna_crossref/` cross-reference set are identical between the two repos.
- Same key limitation: root ground controls n=3 per light condition (low interaction power).

## The pivotal difference — DE model parameterization

| | Manuscript A (main-effects) | Manuscript B (interaction) |
|---|---|---|
| DESeq2 model | `~ light + condition` (additive; light as covariate) | `~ light + condition + light:condition` |
| "Spaceflight effect" reported | **2,132 leaf / 2,582 root** (Flight−Ground, averaged over light) | **527 leaf / 704 root** (condition term = effect at reference light level) |
| Thresholds | padj < 0.05, \|LFC\| ≥ 1 | main padj < 0.05 \|LFC\| ≥ 1; interaction padj < 0.1 \|LFC\| ≥ 0.5 |
| Primary thesis | core microgravity footprint | light **modulates** the spaceflight response (compensatory, r = −0.84) |

The two "spaceflight DEG" counts differ because they are **different estimands**, not a
contradiction: A estimates the light-averaged main effect; B's condition coefficient is
the effect at the reference light level. Manuscript B explicitly discusses comparing its
interaction-model condition coefficients against the additive model (its Results ¶ on
compensatory shrinkage).

## Section / figure overlap

| Content | Manuscript A | Manuscript B | Overlap |
|---|---|---|---|
| QC / PCA, volcano / MA | Figs 1–3 | Supp Figs 1–4 | **High** — same plots |
| GO / KEGG enrichment | Figs 4–5 | Supp Figs 5–6, Table 2 | **High** — same method (A: directional main-effect sets; B: interaction sets) |
| Root cell-type cross-reference | Fig 6 | Fig 1d, 7b, Supp 10 | **High** — same Kajala markers; B extends to radial asymmetry + Cortex_ACT2 OR = 26.9 hotspot |
| Phenylpropanoid / hormone / circadian biology | Results 3.4 + Discussion | Figs 3–5, 7 | **Moderate** — same pathways, reframed as interaction |
| Light × condition interaction | brief (repo `results/light_interaction/`, out of manuscript scope) | **entire thesis** | seed → full paper |
| PhysioSpace stress decoding | **absent** | Figs 8–11, Supp 16–19, Data 11–13 | **Unique to B** |
| Compensatory dynamics, organ hotspots, ISS light-design implications | absent | core | **Unique to B** |

## Assessment & recommendation

The two manuscripts overlap substantially on the **DE / enrichment / cell-type
foundation** (Manuscript A Figs 1–6 ≈ Manuscript B's supplementary DE/enrichment layer),
but diverge completely on **thesis and novel analysis**. Manuscript B's reference list
already frames them as a linked series ("[This study – Analysis 1/2/3]").

**Recommended path:** Manuscript B (PhysioSpace / interaction) is the flagship. Manuscript
A (main-effects) becomes either:
- a **cited companion** (the main-effect / core-footprint reference), or
- **folded into B** as a supplementary main-effects section.

Publishing both as standalone full papers risks self-overlap / dual-publication flags
given the shared Figs 1–6 territory. Whichever route, cross-cite explicitly and make the
model difference (additive vs interaction) the stated reason two analyses exist.

## Open follow-ups

- Decide companion-vs-merge before submitting either.
- If kept separate, add a reciprocal citation in each manuscript and a one-line scope
  statement distinguishing them.
- Fill author/affiliation blocks (placeholders in both drafts).
