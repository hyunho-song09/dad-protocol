# Suggested Reviewers

> Five candidate reviewers selected by the corresponding author (2026-05-08). The DAD team confirms it has no prior personal correspondence with any of the five named below. The list is non-binding on the editors.

---

## Selection criteria

Reviewers should jointly cover three competence axes:

1. **Molecular docking and scoring** — GNINA / Vina / DiffDock / cofolding methodology; pose-evaluation metrics; CNN score calibration.
2. **Bacterial structural bioinformatics** — periplasmic binding proteins, MCPs, Crp/Fnr family; AlphaFold 2 / 3 evaluation; topology prediction (Phobius / DeepTMHMM).
3. **Reproducible workflow engineering for omics** — Snakemake / Nextflow pipelines; container-frozen runtimes; FAIR data and licensing.

A reasonable panel is three reviewers covering all three axes (one each), with two backups.

---

## Candidate slots

| Slot | Required expertise | Candidate name | Affiliation | Email | Justification |
|------|--------------------|----------------|-------------|-------|---------------|
| R1 | Protein structure + docking algorithms | Prof. Chaok Seok | Department of Chemistry, Seoul National University, Republic of Korea | <chaok@snu.ac.kr> | Leading docking-algorithm developer (GalaxyTongDock, GalaxyDock); uniquely qualified to assess GNINA scoring assumptions and the zero-parameter design choice. |
| R2 | Scalable sequence/structure analysis (ColabFold co-author) | Prof. Martin Steinegger | Steinegger Lab, Seoul National University, Korea (formerly Max Planck Institute for Developmental Biology) | (public lab email at https://steineggerlab.com or SNU directory) | Co-author of ColabFold which DAD's Stage 4 invokes; uniquely positioned to assess the AlphaFold-DB / ColabFold integration and reproducibility of the structure-prediction stage. |
| R3 | Reproducible workflow engineering / FAIR | Dr. Björn Grüning | Bioinformatics Group, Department of Computer Science, University of Freiburg, Germany | <bjoern.gruening@gmail.com> | Galaxy / Snakemake ecosystem authority; reviews from the reproducibility-first standpoint that Nature Protocols requires for protocol-class manuscripts. |
| R4 | Membrane-protein topology + AlphaFold validation | Prof. Arne Elofsson | Department of Biochemistry and Biophysics, Stockholm University, Sweden | <arne.elofsson@dbb.su.se> | His group co-developed the Phobius lineage that DAD's optional triage path uses; uniquely qualified to assess the topology-aware pre-docking biological triage. |
| R5 | Metabolomics community / GNPS | Prof. Pieter Dorrestein | Skaggs School of Pharmacy and Pharmaceutical Sciences, UC San Diego, USA | <pdorrestein@ucsd.edu> | Metabolomics community thought leader (GNPS); evaluates whether the many-to-many protein–metabolite matrix output is useful for downstream hypothesis generation in real metabolomics workflows. |

---

## Conflicts of interest to disclose at submission

The corresponding author is asked to flag, per reviewer:

- Any co-authorship within the past 5 years.
- Any current or recent collaborative grant.
- Any current or pending patent overlap.
- Any institutional / departmental affiliation overlap.

A reviewer with any of the above should be moved to "non-preferred" rather than "suggested".

---

## Reviewers to avoid (placeholder)

`[NON_PREFERRED_REVIEWERS_PLACEHOLDER]` — to be filled by the corresponding author. The DAD team does not pre-suggest non-preferred reviewers; this is solely the corresponding author's responsibility.

---

## Notes

- The DAD team has no prior personal correspondence with any of the five candidate reviewers named in this list.
- The suggested-reviewer list is **non-binding** on the editors; Nature Protocols selects reviewers at its discretion.
