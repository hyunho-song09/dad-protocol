# DAD Publication — Placeholder Resolution Guide

> Mr_Prof, 2026-05-08 (Codex Round 8 P0 follow-up).
> One-page operator guide. **22 placeholder hits in `Publication/`** are user-decision items; this document tells the corresponding author exactly what to type, where it goes, and which placeholders are required for which submission mode.

---

## 1. Two submission modes — different placeholder requirements

DAD's evidence base supports two distinct submission entry points. The required-placeholder set is different.

| Submission mode | When to use | Required placeholders | Recommended placeholders |
|-----------------|-------------|----------------------|--------------------------|
| **A. Presubmission enquiry (current default)** | Primary paper not yet accepted *or* author wants editorial confirmation before full upload. Cover letter + manuscript v1 + Tier 1 + RCSB seed evidence is enough. | Categories 2, 3, 4 (GitHub repo, authorship, reviewers) **as placeholders is acceptable** because the submission system does not enforce them. Cover letter must declare conditional status. | Category 1 (primary paper) at minimum should name the paper and its current status (e.g. "in revision at *Journal X*"). |
| **B. Full Nature Protocols submission** | Primary paper accepted or co-submission permitted by editor. | **All 22 placeholders must resolve.** No conditional language remains in cover letter. | — |

This guide assumes mode A unless the user says otherwise.

---

## 2. Six placeholder categories (22 total)

### Category 1 — Supporting primary research paper (4 placeholders)

| Placeholder | What to provide | Where it lands |
|-------------|----------------|----------------|
| `[PRIMARY_PAPER_DOI]` | Published DOI; or `bioRxiv:10.1101/...`; or string `in revision` if not yet accepted | `manuscript/DAD_manuscript_v1.md` Reporting summary §11; `submission_metadata/cover_letter.md`; `submission_metadata/data_availability.md`; `references/references.bib` `PrimaryPaperPlaceholder` entry |
| `[PRIMARY_PAPER_CITATION]` | Full Vancouver-style or APA citation | `manuscript_skeleton.md` §13 References; `cover_letter.md`; `references.bib` |
| `[PRIMARY_PAPER_FIGURE_REF]` | The exact in-paper identifier of the docking table (Reviewer #2 response calls it "Table X") | `manuscript_skeleton.md` §10.1; `06_Report/Mr_Prof/primary_paper_relationship.md` §2 |
| `[PRIMARY_PAPER_SUPPLEMENT_REF]` | "Main text §X" or "Supplementary File N" | `primary_paper_relationship.md` §2 |

**For mode A**: at minimum supply `[PRIMARY_PAPER_DOI] = "in revision"` (or similar) so the cover letter and reporting summary read coherently.

### Category 2 — Code repository (3 placeholders)

| Placeholder | What to provide | Where it lands |
|-------------|----------------|----------------|
| `[GITHUB_REPO_URL_PLACEHOLDER]` | Final public GitHub URL, e.g. `https://github.com/<org>/<repo>` | `code_availability.md`; `data_availability.md`; cover letter |
| `[ZENODO_DOI_PLACEHOLDER]` | Zenodo archive DOI (mint at acceptance, e.g. `10.5281/zenodo.123456`) | `code_availability.md`; `data_availability.md`; reporting summary |
| `[GIT_TAG_PLACEHOLDER]` | Submission-time tag, e.g. `v1.0.0-presubmission` | `code_availability.md` |

**For mode A**: a private repo URL is acceptable as long as the cover letter notes that it will be made public at acceptance.

### Category 3 — Authorship and submission metadata (3 placeholders)

| Placeholder | What to provide | Where it lands |
|-------------|----------------|----------------|
| `[CORRESPONDING_AUTHOR_PLACEHOLDER]` | Full name + email + affiliation of the corresponding author | `cover_letter.md`; `presubmission_enquiry_summary.md`; `conflict_of_interest.md` |
| `[AUTHOR_LIST_PLACEHOLDER]` | Comma-separated full author list with affiliations | `cover_letter.md`; `presubmission_enquiry_summary.md` |
| `[SUBMISSION_DATE_PLACEHOLDER]` | ISO date `YYYY-MM-DD` | `cover_letter.md` |

**Note**: the authorship arrangement with the supporting primary paper's authors must be settled before this category can be filled. See `06_Report/Mr_Prof/primary_paper_relationship.md` §6.

### Category 4 — Reviewer slots (5 placeholders)

| Placeholder | What to provide | Where it lands |
|-------------|----------------|----------------|
| `[REVIEWER_R1_NAME_PLACEHOLDER]` | Docking + CNN scoring expert; full name + affiliation + email | `submission_metadata/suggested_reviewers.md` |
| `[REVIEWER_R2_NAME_PLACEHOLDER]` | Bacterial structural bioinformatics expert | same |
| `[REVIEWER_R3_NAME_PLACEHOLDER]` | Reproducible-pipeline / FAIR expert | same |
| `[REVIEWER_R4_NAME_PLACEHOLDER]` | Backup; topology prediction | same |
| `[REVIEWER_R5_NAME_PLACEHOLDER]` | Backup; metabolomics-coupled inference | same |

The DAD team intentionally does **not** suggest names. Selection is the corresponding author's responsibility, including the COI check listed in `suggested_reviewers.md`.

### Category 5 — Supplementary Information bibliography (5 placeholders)

Five citations are referenced inside SI-Methods §SI-M.1 (triage rules) but not yet in `references.bib`. Add the BibTeX entries before submission.

| Citation | Use site | Suggested key |
|----------|----------|---------------|
| Nielsen, H. *Predicting subcellular localization of proteins by bioinformatic algorithms.* Curr. Top. Microbiol. Immunol. 2017 | SI-M.1 R2 (signal peptide biology) | `Nielsen2017` |
| Drew, D. & Boudker, O. *Shared molecular mechanisms of membrane transporters.* Annu. Rev. Biochem. 2016, 85, 543–572 | SI-M.1 R3 (polytopic IM exclusion) | `DrewBoudker2016` |
| Briegel, A. et al. *Universal architecture of bacterial chemoreceptor arrays.* eLife 2014, 3, e02181 | SI-M.1 R3 path C (MCP architecture) | `Briegel2014` |
| Zhulin, I. B. *The superfamily of chemotaxis transducers.* Adv. Microb. Physiol. 2001, 45, 157–198 | SI-M.1 R5 (sensor evolution) | `Zhulin2001` |
| Krell, T. et al. *Bacterial sensor kinases.* Annu. Rev. Microbiol. 2010, 64, 539–559 | SI-M.1 R5 (sensor / SBP boost rationale) | `Krell2010` |

**For mode A**: these can stay as placeholder citations as long as the SI-Methods file uses them as `[Nielsen 2017]` text-citation form. **For mode B**: BibTeX entries must be present.

### Category 6 — Embargo / co-submission (2 placeholders, decision items)

| Decision | What to record | Where it lands |
|----------|---------------|----------------|
| Embargo timing | "Hold until primary paper acceptance"; "Co-submit on Y date"; or "No embargo" | Cover letter §3 |
| Co-submission consent | Editor's response confirming acceptance of co-submission status | Internal record; reference only in cover letter if affirmative |

These are decision items rather than typed strings; once decided, the cover letter is updated accordingly.

---

## 3. One-shot user input form

Filling the form below in a single user turn provides everything needed to convert the package from mode A to mode B.

```yaml
# DAD Publication — placeholder resolution
# Fill each field; leave as ~null~ to keep placeholder for mode A presubmission.

primary_paper:
  doi:                  # e.g. "10.1038/s41xxx-xxx"  or  "in revision"
  citation:             # full citation string
  figure_ref:           # e.g. "Table 4"  or  "Supplementary Table S3"
  supplement_ref:       # e.g. "Main text §3.2"  or  "Supplementary File 2"

repository:
  github_url:           # e.g. "https://github.com/<org>/dad"
  zenodo_doi:           # e.g. "10.5281/zenodo.NNNNNN"  (mint at acceptance)
  git_tag:              # e.g. "v1.0.0-presubmission"

authorship:
  corresponding_author: # "Name <email@institution>"
  author_list:          # "Surname, F. Initial; ... ; Last, F. Initial"
  submission_date:      # YYYY-MM-DD

reviewers:
  R1_name:              # docking + CNN scoring expert
  R2_name:              # bacterial structural bioinformatics expert
  R3_name:              # reproducible pipeline expert
  R4_name:              # backup — topology prediction
  R5_name:              # backup — metabolomics inference

si_bibliography:
  nielsen_2017:         # full citation
  drew_boudker_2016:    # full citation
  briegel_2014:         # full citation
  zhulin_2001:          # full citation
  krell_2010:           # full citation

decisions:
  embargo:              # "hold-until-acceptance" | "co-submit:YYYY-MM-DD" | "none"
  co_submission_status: # "requested" | "approved" | "declined" | "n/a"
```

Once the user fills this form, Mr_Prof can substitute every placeholder across the package in a single round.

---

## 4. Quick verification

Two grep commands check resolution status at any time:

```bash
# All remaining placeholders in Publication/
grep -RE "\[[A-Z_]+_PLACEHOLDER\]|\[PRIMARY_PAPER_[A-Z_]+\]|\[GITHUB_REPO_URL\]|\[ZENODO_DOI\]|\[GIT_TAG\]" Publication/

# Mode-A-blocking only (Category 1):
grep -RE "\[PRIMARY_PAPER_(DOI|CITATION|FIGURE_REF|SUPPLEMENT_REF)\]" Publication/
```

Mode A is satisfied when the second command's output is empty (or contains only `"in revision"`-style explicit-conditional strings). Mode B is satisfied when the first command's output is empty.

---

*This guide is the canonical reference for the corresponding author. Internal categorisation cross-references `06_Report/Mr_Prof/manuscript_skeleton.md` §15 (open questions) and `06_Report/Mr_Prof/primary_paper_relationship.md` §7 (primary paper open items).*
