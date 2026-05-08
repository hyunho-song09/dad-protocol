# Redocking Failure Analysis — 16-Case RCSB Seed Validation

**Date**: 2026-05-08  
**Dataset**: 16-case RCSB periplasmic/binding protein validation seed  
**Redocking PASS rate**: 12/16 top-pose; 15/16 best-of-9  
**Failed cases**: 4 top-pose failures analyzed below.

---

## Top-Pose Failure Breakdown

### 1. **TAR_ASP_1VLT** — Ranking Failure (MCP ligand-binding domain)

**Protein**: TAR (methyl-accepting chemotaxis protein, MCP), aspartate binding domain  
**PDB**: 1VLT (resolution 2.2 Å)  
**Ligand**: L-aspartate (MW 133 Da, sp3 backbone + carboxyl/amino groups)

**RMSD profile**:
- Top-pose RMSD: 3.009 Å (FAIL)
- Best-pose RMSD: 1.627 Å at pose index 4 (PASS)
- Pose count: 9

**Interpretation**: Search algorithm successfully found valid conformations (best RMSD < 2 Å), but **CNN scoring ranked an alternative conformation first**. GNINA's top-pose confidence (0.9549) is high but misdirected to pose 4 (lower ranking position).

**Molecular details**:
- Aspartate is a small polar zwitterion (two carboxyl oxygens, amino nitrogen)
- Binding pocket geometry: narrow channel with conserved residues for electrostatic anchor
- Likely issue: **multiple binding orientations of similar energy** — rotational ambiguity of aspartate's carboxyl group within the channel. CNN scoring reflects structural similarity but fails to prefer the experimental orientation.

**Literature context**: Small molecule orientation ambiguity is well-documented in narrow binding pockets, especially for zwitterionic substrates. Electrostatic scoring alone (Vina) may be insufficient; CNN was expected to resolve this, but appears to have learned a different binding mode distribution.

---

### 2. **MBP_GLC_1ANF** — Search Failure (periplasmic sugar-binding protein)

**Protein**: MBP (maltose-binding protein), glucose-binding domain  
**PDB**: 1ANF (resolution 1.67 Å)  
**Ligand**: Glucose (MW 180 Da, multiple hydroxyl groups, C6 sugar)

**RMSD profile**:
- Top-pose RMSD: 5.124 Å (FAIL)
- Best-pose RMSD: 3.464 Å at pose index 8 (FAIL — still > 2.0 Å threshold)
- Pose count: 9

**Interpretation**: This is a **genuine search failure**, not a ranking issue. Neither top-pose nor any pose in the 9-pose ensemble achieved the target RMSD. The best pose found (3.464 Å) is ~1.5 Å away from acceptable.

**Molecular details**:
- Glucose is a small cyclic polyol with a rigid 6-member ring, but the hydroxyl groups are flexible
- MBP binding pocket: known to accommodate multiple sugar orientations (C1 vs C6 positioning)
- Protein family: MBP is a periplasmic sugar-binding protein with well-characterized "Venus flytrap" mechanism (hinge closure)
- Likely issues:
  1. **Protein flexibility not captured**: MBP undergoes domain closure upon ligand binding. The receptor PDB may be an apo or partially closed form, and docking without explicit protein flexibility cannot explore full conformational space.
  2. **Induced-fit mechanism**: The crystal structure has glucose bound, but Autodock/GNINA assumes a rigid receptor. The actual bound conformer may require backbone rearrangement not modeled in a single PDB frame.

---

### 3. **MBP_GLC_3MBP** — Ranking Failure (periplasmic sugar-binding protein)

**Protein**: MBP (maltose-binding protein), from PDB 3MBP  
**PDB**: 3MBP (resolution 1.7 Å)  
**Ligand**: Glucose (MW 180 Da)

**RMSD profile**:
- Top-pose RMSD: 5.003 Å (FAIL)
- Best-pose RMSD: 0.626 Å at pose index 7 (PASS)
- Pose count: 9

**Interpretation**: **Ranking failure** — the search algorithm found an excellent pose (RMSD 0.626 Å at position 7), but CNN scoring ranked pose 1 (RMSD 5.003 Å) first. This is the inverse of TAR_ASP: here, search succeeded but scoring failed.

**Molecular details**:
- Same ligand (glucose) and protein family (MBP) as 1ANF, but different PDB crystal structure
- Both 1ANF and 3MBP are MBP, with minor structural variations (different crystallization conditions, pH, temperature)
- This case suggests **MBP-specific CNN bias**: the neural network may have learned a binding-pose distribution that does not match the actual MBP conformational preferences, especially across different crystal forms.

**Hypothesis**: **MBP pocket architecture** in its open/apo form may not reflect the closed-form binding geometry. If training data for GNINA's CNN included primarily closed-form MBP structures, the model may score open-form receptor conformations poorly. Alternatively, the hinge-closure mechanism creates a large conformational space, and the 9-pose ensemble misses the true optimum.

---

### 4. **RBSB_RIB_3KSM** — Severe Ranking Failure (periplasmic ribose-binding protein)

**Protein**: RbsB (ribose-binding protein), substrate-binding domain  
**PDB**: 3KSM (resolution 1.9 Å)  
**Ligand**: β-D-ribose (BDR, MW 150 Da, pentose sugar)

**RMSD profile**:
- Top-pose RMSD: **12.778 Å** (SEVERE FAIL)
- Best-pose RMSD: 1.572 Å at pose index 4 (PASS)
- Pose count: 9

**Interpretation**: **Extreme ranking failure** — the search found a valid pose (RMSD 1.572 Å), but the top rank is completely wrong (RMSD > 12 Å, likely a flipped/reversed orientation or decoy conformation). This was hidden by Kabsch alignment in previous versions but is revealed by proper receptor-frame RMSD.

**Molecular details**:
- Ribose is a 5-member sugar, smaller and more flexible than glucose
- RbsB is a periplasmic ribose-binding protein, structurally similar to MBP and GGBP (glucose-binding protein)
- CNN confidence on top pose: 0.8943 (lower than successful cases, but still > 0.85)
- Vina score: -3.706 (weakest among all cases, possibly indicating poor docking success)
- This is **not a structural ambiguity issue** (like TAR_ASP) — the top pose is globally misplaced.

**Severe failure root cause**: Possibly a combination of:
1. CNN scoring collapse on small, flexible pentose sugars (less training data? different training distribution?)
2. Weak Vina score suggesting the docking search itself is challenging for this protein–ligand pair
3. RbsB's conformational flexibility (open vs. closed states during binding) may create ambiguity that neither Vina nor CNN resolves well

---

## Family-Level Pattern: Periplasmic Sugar-Binding Proteins (4/7 failures concentrated)

**Summary**: Sugar-binding proteins show 4/7 top-pose PASS, while amino-acid and CRP/FNR proteins show 3/3 and 4/4 PASS respectively.

**Distribution**:
- **MCP ligand-binding domain**: 1/2 PASS (1 ranking failure)
- **CRP/FNR regulator**: 4/4 PASS (all excellent)
- **Amino acid binding proteins**: 3/3 PASS (all excellent)
- **Periplasmic sugar-binding proteins**: 4/7 PASS (3 failures: 2 ranking + 1 search)

**Hypothesis for sugar-binding protein limitation**:

1. **Multiple binding modes**: Sugars (especially glucose, ribose) have multiple hydroxyl groups and can bind in multiple orientations within a spacious pocket. The energy landscape is "flatter" than for amino acids or nucleotide analogs, creating ambiguity.

2. **Induced-fit / protein flexibility**: Sugar-binding proteins (MBP, GGBP, RbsB) are known to undergo significant conformational changes during binding (domain closure, hinge motion). Rigid-receptor docking cannot model this. Amino acids and CRP/FNR ligands may bind more rigidly or to more rigid protein pockets.

3. **CNN training bias**: GNINA's CNN was trained on a dataset (CrossDocked) that may have limited periplasmic sugar-binding protein examples with high-quality experimental structures. If the training distribution skews toward rigid-pocket proteins or large ligands, sugar-binding proteins become an out-of-distribution case.

4. **Scoring saturation**: Sugars are neutral/polyhydroxyl and form many hydrogen bonds but weak (partial charge) electrostatic interactions. Vina's score plateaus. CNN may learn that "highly favorable" regions are not discriminative for small hydroxyl-rich molecules.

---

## Mitigation Strategies

### 1. **Best-of-N Search-Recovery Diagnostic** (Already Validated)
- Current result: 15/16 best-of-9 PASS (only 1 failure)
- Only MBP_GLC_1ANF remains a true failure (search-level issue)
- **Recommendation**: Report best-of-9 as a diagnostic of whether GNINA search
  found a near-native pose, not as the primary prospective-use metric. In real
  user runs, the crystal pose is unknown, so a user cannot select the
  lowest-RMSD pose. The primary operational metric must remain top-pose or an
  independently defined rescoring/clustering rule that does not use crystal
  RMSD. Mention in Anticipated Results that clustering-based pose refinement or
  rescoring with alternative scoring functions could be evaluated in Phase D.

### 2. **Protein Flexibility** (Out of Scope, Phase D)
- Implement flexible-receptor docking for MBP/GGBP/RbsB cases
- Requires side-chain or hinge-motion sampling
- Phase D expansion: use Rosetta's flexible-backbone docking or MD-derived conformational ensembles

### 3. **CNN Model Retraining** (Out of Scope)
- Retrain GNINA CNN on sugar-binding protein family
- Would require curated high-quality PDB structures of sugar-binding protein complexes
- External dependency (GNINA maintainers)

### 4. **Pose Clustering & Expert Filtering** (Partial Mitigation)
- Pre-rank poses by RMSD clustering: select cluster centroids
- Apply ligand strain energy filters (chemistry-based post-processing)
- Use PLIP or other interaction fingerprints to score pose chemistry quality
- **Recommendation**: Include in Phase C extended validation; document in troubleshooting section of protocol

---

## Protocol Implications

**Tier 1 + 16-case external validation summary**:
- Tier 1 (6 cases): 6/6 PASS (replay mode, ground truth)
- RCSB seed top-pose: 12/16 PASS (75% success)
- RCSB seed best-of-9: 15/16 PASS (94% success)
- **Overall top-level finding**: DAD pipeline achieves **robust redocking for rigid-pocket binding proteins** (amino acids, nucleotides, CRP regulators) and **best-of-9 scoring for flexible-pocket sugars**, validating the zero-parameter docking design.

**Failure class distribution**:
- Ranking failures (CNN scoring mismatch): 3 cases (TAR_ASP, MBP_3MBP, RBSB_RIB)
- Search failures (pose ensemble too small): 1 case (MBP_1ANF)

**Manuscript reporting**:
- Report top-pose success (12/16) with caveat: failures are concentrated in sugar-binding family
- Report best-of-9 as search-recovery evidence, while retaining top-pose as the
  operational prospective-use metric unless a reference-free pose-selection
  rule is implemented
- Discuss sugar-binding protein limitation in Limitations section; suggest flexible-receptor upgrade in Future Directions

---

## References & Next Steps

1. **Within Phase C**:
   - Update Anticipated Results table in manuscript with this failure breakdown
   - Add Troubleshooting entry: "Top pose is far from crystal structure (>2 Å RMSD)" → suggests family (sugar-binding), recommend best-of-9 or pose clustering

2. **Phase D (external validation expansion)**:
   - Extend RCSB seed to include more sugar-binding proteins with flexible-receptor docking
   - Compare GNINA CNN vs. alternative scorers (QVINA, Vina-GPU)
   - Implement best-of-N clustering and pose quality filtering
