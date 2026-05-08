# Citation Index — DAD Manuscript v1

> Mr_Prof, 2026-05-08. Per-section citation map. Cross-checks `references.bib` against in-text use sites in `Publication/manuscript/DAD_manuscript_v1.md` and `Publication/supplementary/SI_text.md`.

---

## Section → BibTeX key map

| Manuscript section | BibTeX key | Use site (paraphrased) |
|--------------------|-----------|------------------------|
| Introduction §1.2 | Corso2023 | DiffDock single-target assumption |
| Introduction §1.2 | Trott2010 | AutoDock Vina single-target assumption |
| Introduction §1.2 | McNutt2025 | GNINA standalone |
| Introduction §1.2 | Krishna2024RFAA | RFAA cofolding |
| Introduction §1.2 | Wohlwend2024Boltz | Boltz-1 cofolding |
| Introduction §1.2 | Abramson2024 | AlphaFold 3 cofolding |
| Introduction §1.2 | Stark2022EquiBind | EquiBind regression |
| Introduction §1.2 | Hekkelman2022AlphaFill | AlphaFill ligand transplant |
| Introduction §1.3 | Mirdita2022 | ColabFold |
| Introduction §1.3 | Krivak2018 | P2Rank |
| Introduction §1.3 | McNutt2025 | GNINA 1.3 |
| Introduction §1.3 | RDKit | RDKit |
| Introduction §1.3 | Adasme2021PLIP | PLIP |
| Introduction §1.3 | Koster2012Snakemake | Snakemake |
| Development of the protocol §2 | Sung2026KeratinStarvationFI | supporting primary paper (iScience 2026, peer-reviewed published) |
| Development of the protocol §2 | Mirdita2022, McNutt2025 | iteration-2 notebooks |
| Development of the protocol §2 | Hallgren2022DeepTMHMM | DeepTMHMM default |
| Development of the protocol §2 | Kall2007Phobius | Phobius opt-in |
| Comparison §5 (Table 2) | Corso2023, Trott2010, McNutt2025, Hekkelman2022AlphaFill, Stark2022EquiBind, Krishna2024RFAA, Wohlwend2024Boltz | 8-tool table rows |
| Comparison §5 (license-axis text) | Berman2000PDB, Liu2007BindingDB, Zhang2024BioLiP, Francoeur2020, Mendez2019, Varadi2022AFDB, Wishart2022HMDB, Haug2020MetaboLights | reference DB licenses |
| Experimental design §5.1 | Sung2026KeratinStarvationFI | Tier 1 ground-truth source |
| Experimental design §5.4 | Buttenschoen2024 | PoseBusters as PB-validity reference; Round-6 quantitative analysis context |
| Experimental design §5.4 | Francoeur2020 | CrossDocked2020 = GNINA training distribution |
| Experimental design §5.5 | (no external citation; runtime_freeze.md) | — |
| Materials §6.4 | Hallgren2022DeepTMHMM | DeepTMHMM install |
| Materials §6.4 | Kall2007Phobius | Phobius opt-in |
| Procedure §7 (Stage 4) | Mirdita2022, Jumper2021, Abramson2024, Wohlwend2024Boltz | structure-prediction options |
| Procedure §7 (Stage 5) | Krivak2018 | P2Rank |
| Procedure §7 (Stage 8) | McNutt2025 | GNINA 1.3 default |
| Procedure §7 (Stage 9) | Adasme2021PLIP, Hamelryck2003BioPDB | interaction profiling |
| Procedure §7 (Stage 11) | Pettersen2021ChimeraX | ChimeraX `.cxc` |
| Anticipated Results §9.1 | Sung2026KeratinStarvationFI | Tier 1 replay anchor |
| Anticipated Results §9.2 | Berman2000PDB | RCSB seed source |
| Anticipated Results §9.2 | McNutt2025, Francoeur2020 | GNINA 1.3 + CrossDocked2020 training |
| Limitations §10.1 | McNutt2025, Francoeur2020 | within-case CNN ranking limitation |
| Limitations §10.2 | Francoeur2020 | CrossDocked2020 training-distribution bias |
| Limitations §10.3 | Buttenschoen2024 | PoseBusters validity |
| Reporting summary §11 | (all software refs above) | — |

## Supplementary section → BibTeX key map

| SI section | BibTeX key | Use site |
|------------|-----------|----------|
| SI-M.1 R1 length | Jumper2021 | pLDDT meaningful only above ~ 50 aa |
| SI-M.1 R2 SP clip | Kall2007Phobius | signal peptide / TM topology backbone |
| SI-M.1 R3 (DeepTMHMM default) | Hallgren2022DeepTMHMM | nTM topology classification |
| SI-M.1 R3 nTM ≥ 7 | DrewBoudker2016Transporters | polytopic IM exclusion + induced-fit context |
| SI-M.1 R3 path C | Briegel2014ChemoreceptorArchitecture | MCP cytoplasmic-domain architecture |
| SI-M.1 R4 | Krivak2018 | P2Rank length sensitivity |
| SI-M.1 R5 | Krell2011Diversity | chemoreceptor functional-class boost rationale |
| SI-M.3 RCSB seed | Berman2000PDB | source database |
| SI-M.4 RMSD evaluator | RDKit | symmetry-corrected CalcRMS |
| SI-M.5 score parser | McNutt2025 | GNINA SDF property tags |
| SI-Notes (failure analysis) | (no new citations; cross-references main text) | — |
| SI-License Table SI-L | Berman2000PDB, Varadi2022AFDB, Liu2007BindingDB, Zhang2024BioLiP, Francoeur2020, Mendez2019, Haug2020MetaboLights, Wishart2022HMDB, Buttenschoen2024 | reference DB licenses |

## Coverage check

- All 33 BibTeX entries currently in `references.bib` are cited at least once above.
- The five SI-Methods placeholders are now resolved with team-lead-selected citations: `Kall2007Phobius`, `Hallgren2022DeepTMHMM`, `Krell2011Diversity`, `DrewBoudker2016Transporters`, `Briegel2014ChemoreceptorArchitecture`. (User-input form 2026-05-08.)
- `Sung2026KeratinStarvationFI` (iScience 2026, Cell Press; PII S2589-0042(26)01130-2) replaces `PrimaryPaperPlaceholder`. Status: peer-reviewed and published; Nature Protocols' supporting-research-paper requirement satisfied.

## Reference-count summary

| Body section | Refs cited |
|--------------|-----------|
| Introduction | 13 |
| Development of the protocol | 5 |
| Applications | 0 (cross-references) |
| Comparison | 8 (table rows) + 8 (license-axis) |
| Experimental design | 3 |
| Materials | 2 |
| Procedure | 9 |
| Troubleshooting | 0 (cross-references) |
| Anticipated Results | 4 |
| Limitations | 3 |
| Reporting summary | n/a (lists software) |
| **Body unique total** | ~ 22 (after deduplication) |
| **+ SI-Methods unique** | + 5 (Kall2007Phobius, Hallgren2022DeepTMHMM, Krell2011Diversity, DrewBoudker2016Transporters, Briegel2014ChemoreceptorArchitecture) |
| **+ Supporting primary** | + 1 (Sung2026KeratinStarvationFI) |
| **Total in `references.bib`** | **33 entries** |
| **Submission target** | 25–30 references (Nature Protocols typical range) — final pruning at submission |
