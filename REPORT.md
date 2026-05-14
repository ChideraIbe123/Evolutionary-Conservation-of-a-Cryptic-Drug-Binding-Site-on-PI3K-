# Evolutionary Conservation of a Cryptic Drug Binding Site on PI3Kα

**Chidera Ibe** — coibe2@illinois.edu

---

## Abstract

Phosphatidylinositol 3-kinase alpha (PI3Kα) is mutated in roughly 40% of hormone receptor-positive breast cancers, making it one of the most important oncology drug targets. In 2024, D.E. Shaw Research reported RLY-2608, a first-in-class inhibitor that binds a cryptic allosteric pocket on PI3Kα and selectively inhibits the mutant form while sparing wild-type signaling in healthy tissue. This project asks whether the pocket residues are evolutionarily conserved, and whether conservation can predict which residues develop resistance mutations in patients. Using from-scratch implementations of Needleman-Wunsch with affine gap penalties, progressive multiple sequence alignment, and the neighbor-joining phylogeny algorithm, I aligned 20 vertebrate PI3Kα orthologs and the four human PI3K paralogs (α, β, γ, δ), then mapped per-residue Shannon entropy conservation scores onto the co-crystal structure PDB 8TSD. Three findings stand out. First, PI3Kα is under extreme purifying selection — 87% of its 1068 residues are perfectly conserved across vertebrates. Second, the cryptic pocket is even more conserved than the protein average (mean conservation 0.985 vs. 0.969). Third, every clinically observed resistance mutation from Varkaris et al. (2024) occurs at a position that is 100% conserved across the species sampled. The corollary is that conservation, by itself, does not discriminate resistance-prone positions from the rest of the protein, because most of the protein is maximally conserved. But the deeper finding is that resistance at these hyper-conserved sites means tumors are paying a functional fitness cost to escape the drug — a clinically encouraging observation that reframes how we think about resistance in allosteric pockets of highly constrained enzymes.

## 1. Introduction

PI3Kα, encoded by the *PIK3CA* gene, is a lipid kinase that phosphorylates phosphatidylinositol 4,5-bisphosphate to generate PIP3, activating AKT/mTOR signaling that drives cell growth and survival. Gain-of-function mutations in *PIK3CA* — particularly H1047R, E545K, and E542K — are among the most common oncogenic drivers in breast cancer. First-generation PI3K inhibitors such as alpelisib were approved for PIK3CA-mutant breast cancer but cause serious on-target toxicity (hyperglycemia, rash) because they inhibit wild-type PI3Kα in healthy cells. RLY-2608 addresses this by binding a cryptic allosteric pocket at the C2–kinase domain interface that preferentially forms in mutant PI3Kα. The drug is selective for the alpha isoform and favors the mutant over wild-type, and is now in clinical trials. In 2024, Varkaris et al. (*Cancer Discovery*) reported secondary resistance mutations observed in patients treated with PI3Kα inhibitors, including W780R, Q859K/H, E726K, and I817F — positions outside the orthosteric ATP site but clustered around the allosteric pocket. Whether these resistance positions differ in evolutionary conservation from the rest of the protein has not, to my knowledge, been systematically analyzed. If the pocket is highly conserved, resistance is constrained because mutations there likely come at a fitness cost. If the pocket is variable, conservation might predict which positions are most vulnerable. This project tests both hypotheses by computing per-residue conservation across vertebrate PI3Kα orthologs and mapping the scores onto the 8TSD co-crystal structure.

## 2. Materials and Methods

### 2.1 Datasets

Four datasets were used. The co-crystal structure PDB 8TSD (2.70 Å resolution) provides the 3D coordinates of PI3Kα bound to RLY-2608. The ligand residue name is XUZ. PI3Kα ortholog sequences were retrieved from the Ensembl REST API via a homology query on human *PIK3CA* (UniProt P42336, 1068 aa); of the 222 orthologs returned, I selected 20 representative vertebrate species spanning primates, rodents, carnivores, ungulates, marsupials, birds, reptiles, amphibians, and fish. Paralog sequences for the four human PI3K catalytic isoforms (PIK3CA P42336, PIK3CB P42338, PIK3CG P48736, PIK3CD O00329) were fetched from UniProt. Resistance mutation positions were taken from Varkaris et al., *Cancer Discovery* 2024.

### 2.2 Sequence Alignment (Notebook 1)

Needleman-Wunsch with affine gap penalties was implemented from scratch in NumPy. The algorithm uses three dynamic programming matrices — M, Ix, Iy — corresponding to match/mismatch, gap-in-y, and gap-in-x states, with gap-open penalty -10 and gap-extend penalty -0.5 on the BLOSUM62 substitution matrix. Traceback includes boundary guards to handle the case when one index reaches zero. Correctness was verified by comparing alignment scores against Biopython's `PairwiseAligner` on the same sequences; the scores agreed within floating-point tolerance. Progressive multiple sequence alignment was implemented using profile-profile alignment following a guide tree (neighbor-joining, see §2.3). At each internal node of the tree, two child profiles are merged using a modified Needleman-Wunsch where column scores are computed as frequency-weighted averages of all pairwise BLOSUM62 scores between residues in the two profile columns. The final ortholog MSA contains 20 sequences over 1073 columns; the paralog MSA contains 4 sequences over 1206 columns.

### 2.3 Phylogenetic Tree (Notebook 2)

An all-vs-all pairwise distance matrix was built from Needleman-Wunsch alignment percent identities, converted to evolutionary distances with the Jukes-Cantor correction for proteins. The neighbor-joining algorithm of Saitou and Nei (1987) was implemented from scratch: (i) compute net divergence for each taxon, (ii) compute the Q-matrix, (iii) find the minimum Q entry and join that pair into a new node, (iv) compute branch lengths, (v) update the distance matrix with distances to the new node, (vi) repeat until two taxa remain. The resulting Newick tree was saved and used as the guide tree for progressive MSA in Notebook 1. Topology was verified by confirming expected clades: human-chimp are nearest neighbors, mouse-rat are nearest neighbors, chicken-turkey cluster together, fish form an outgroup to tetrapods.

### 2.4 Conservation Scoring (Notebook 3)

Per-column Shannon entropy was computed from the ortholog MSA as H = -Σ p_i log₂(p_i), where p_i is the frequency of amino acid i in the column (gaps excluded). Entropy was normalized to a conservation score as 1 - H / log₂(20), so that a score of 1.0 indicates perfect conservation and 0 indicates maximum entropy. Columns were mapped back to human PI3Kα residue numbering by indexing only non-gap positions of the human reference row.

### 2.5 Pocket Identification and Structural Mapping

PDB 8TSD was parsed with Biopython's PDBParser, restricting analysis to chain A (the p110α catalytic subunit). The RLY-2608 ligand was identified by its HETATM residue name (XUZ). Using Biopython's NeighborSearch, all chain A protein residues with any atom within 5.0 Å of any ligand atom were defined as the cryptic pocket. Conservation scores were written into the B-factor column of a modified PDB file and visualized in Jupyter with py3Dmol, rendering the cartoon colored by conservation (blue = conserved, red = variable), pocket residues as sticks, and the ligand in yellow.

### 2.6 Statistical Analysis (Notebook 4)

Three groups of residue positions were compared: (i) resistance mutation sites from Varkaris et al. (W780, Q859, E726, I817), (ii) primary oncogenic sites (H1047, E545, E542), and (iii) all other residues. Distributions of conservation scores across groups were compared with the Mann-Whitney U test. A permutation test (10,000 iterations) was run by sampling random residue sets of size equal to the resistance group and computing the null distribution of mean conservation. Fisher's exact test was applied to a 2×2 contingency table of resistance × pocket membership. All statistical tests were performed in SciPy.

## 3. Results

### 3.1 PI3Kα is under extreme purifying selection

Across 20 vertebrate species spanning approximately 450 million years of evolution (human through coelacanth and ray-finned fish), 87.3% of the 1068 residues in human PI3Kα have a conservation score greater than 0.9. The median conservation score is 1.0, the mean is 0.969, and the minimum is substantially above zero. This is consistent with PI3Kα being a non-redundant signaling hub whose function cannot tolerate most substitutions.

### 3.2 The cryptic pocket is more conserved than the protein average

The NeighborSearch analysis identified 26 residues within 5 Å of the RLY-2608 ligand in PDB 8TSD, spanning positions 808–1025. Mean conservation at pocket residues is 0.985, compared to 0.969 for the rest of the protein. The pocket is not a randomly exposed crevice — its residues appear to be structurally constrained, likely because they mediate C2–kinase domain contacts that are critical for kinase regulation even when the allosteric pocket is closed.

### 3.3 Resistance mutations occur at perfectly conserved positions

All four clinically observed resistance mutation positions — W780, Q859, E726, and I817 — score 1.0 for conservation (every amino acid in the 20-species alignment is identical at these columns). The same is true for the three primary oncogenic positions (H1047, E545, E542). This means cancer-associated mutations arise at sites where evolution has not tolerated *any* variation across vertebrate history.

### 3.4 Conservation alone does not discriminate resistance-prone sites

Because both the resistance positions and the bulk of the protein are essentially saturated at conservation 1.0, the Mann-Whitney U test comparing resistance vs. non-resistance residues returns p = 0.28, and the permutation test p-values are similarly non-significant (p = 0.36 for "resistance more conserved than random"). The biological interpretation is not that conservation is irrelevant, but that conservation has insufficient discriminating power in a protein where nearly every position is conserved. Supplementary information — structural context, proximity to the binding pocket, contact with the ligand, participation in domain-domain interfaces — is needed to refine predictions.

### 3.5 Isoform divergence explains drug selectivity

While the four human PI3K catalytic isoforms (α, β, γ, δ) share the same overall domain architecture, their pairwise sequence identities are only 29–56%: PI3Kα vs PI3Kβ 38.6%, PI3Kα vs PI3Kδ 38.9%, PI3Kβ vs PI3Kδ 56.4%, and PI3Kγ is most divergent from all others at ~29%. Of the 26 pocket-lining residues in PI3Kα, only 15 are identical across all four isoforms; 11 differ in at least one other paralog. This divergence is sufficient to explain why RLY-2608 inhibits PI3Kα selectively over β/γ/δ — the pocket has enough α-specific residues to serve as a recognition signature, even though the pocket is conserved across vertebrate α orthologs.

## 4. Discussion

The headline finding is that the RLY-2608 cryptic pocket sits in a region of PI3Kα under extreme evolutionary constraint, and that clinical resistance mutations are not occurring at variable positions that evolution permitted to drift, but at perfectly conserved positions where evolution has refused to accept any change. This inverts the naive hypothesis that resistance should track conservation. The correct interpretation is that strong drug-imposed selective pressure can overcome even deep evolutionary constraint, pushing tumors to accept mutations that are otherwise forbidden.

There are three reasons this is clinically encouraging rather than discouraging. First, mutations at conserved sites likely carry a functional fitness cost — a tumor that acquires W780R or Q859K to escape RLY-2608 may be sacrificing kinase activity, regulatory contacts, or stability, which could render it vulnerable to other interventions or impair its proliferative advantage. Second, because the pocket is so conserved, the number of positions that can even structurally accommodate a resistance mutation is limited — the mutational landscape available to the tumor is narrow. Third, the pocket is simultaneously conserved across orthologs and divergent across paralogs, which means the pharmacology — selectivity plus constrained resistance — is a feature of the binding site design, not a coincidence.

There are also limitations. First, the 20-species ortholog set is all vertebrate, and for a gene as conserved as PI3Kα, vertebrate variation may be insufficient to reveal finer-grained differences that might distinguish resistance-prone from resistance-resistant positions. Including invertebrate and fungal homologs (where alignment becomes harder but evolutionary signal is richer) could provide more discriminating conservation estimates. Second, the resistance mutation list from Varkaris et al. is small (five mutations at four positions), limiting statistical power. As more patients are treated with PI3Kα inhibitors, the resistance mutation catalog will grow and the analysis can be repeated. Third, conservation is one of several features that might predict resistance; combining conservation with structural features (solvent accessibility, distance to ligand, role in allosteric activation) in a multivariate model would likely outperform any single feature.

A methodological observation is worth noting. Implementing Needleman-Wunsch, progressive MSA, and neighbor-joining from scratch confirmed that these textbook algorithms produce biologically reasonable outputs — the neighbor-joining tree correctly recovered vertebrate phylogeny (primates clade together, rodents clade together, birds clade together, fish form an outgroup), and the progressive MSA sum-of-pairs scores were consistent with a high-quality alignment. The boundary conditions in the affine-gap Needleman-Wunsch traceback are not entirely obvious; a careful implementation must handle the cases where one index reaches zero before the other, which is easy to overlook in a naive port of the standard pseudocode.

## 5. Conclusion

Using only algorithms from a bioinformatics course curriculum — Needleman-Wunsch with affine gaps, progressive MSA guided by a neighbor-joining tree, and Shannon-entropy conservation scoring — this project produced a structurally and clinically interpretable analysis of the RLY-2608 cryptic drug pocket on PI3Kα. The pocket is more conserved than the protein average; every known clinical resistance mutation site is perfectly conserved across vertebrates; and divergence between the four PI3K paralogs at the pocket explains the drug's isoform selectivity. Conservation cannot, by itself, predict resistance because PI3Kα is too conserved globally — but the fact that resistance mutations occur at fully conserved sites is itself informative, indicating that tumors escaping RLY-2608 are paying a functional cost, which is a clinically favorable property of the drug.

## References

- Saitou N, Nei M (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution* 4(4):406–425.
- Needleman SB, Wunsch CD (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. *Journal of Molecular Biology* 48(3):443–453.
- Varkaris A, et al. (2024). Discovery and clinical proof-of-concept of RLY-2608, a first-in-class mutant-selective allosteric PI3Kα inhibitor that decouples anti-tumor activity from hyperinsulinemia. *Cancer Discovery* 14(2):240–257.
- Henikoff S, Henikoff JG (1992). Amino acid substitution matrices from protein blocks. *Proceedings of the National Academy of Sciences* 89(22):10915–10919.
- PDB entry 8TSD — Co-crystal structure of PI3Kα bound to RLY-2608. RCSB Protein Data Bank.
- UniProt entries: P42336 (PIK3CA), P42338 (PIK3CB), P48736 (PIK3CG), O00329 (PIK3CD).
- Ensembl Rest API — https://rest.ensembl.org

## Appendix: Repository Structure

```
Bio Project/
├── requirements.txt                 # Python dependencies
├── src/
│   ├── fetchers.py                 # Ensembl / UniProt / RCSB data fetching
│   ├── utils.py                    # BLOSUM62 loader, FASTA I/O
│   ├── alignment.py                # Needleman-Wunsch, progressive MSA, UPGMA
│   ├── tree.py                     # Neighbor-joining (from scratch)
│   └── conservation.py             # Shannon entropy scoring
├── notebooks/
│   ├── 01_alignment.ipynb
│   ├── 02_phylogenetic_tree.ipynb
│   ├── 03_conservation_structure.ipynb
│   └── 04_resistance_analysis.ipynb
├── data/
│   ├── raw/                        # Cached API downloads
│   ├── processed/
│   │   ├── alignments/             # Ortholog & paralog MSAs (FASTA)
│   │   └── trees/                  # Neighbor-joining tree (Newick)
│   └── results/
│       ├── conservation_scores.csv
│       ├── pocket_residues.csv
│       ├── final_analysis.csv
│       ├── phylogenetic_tree.png
│       ├── conservation_landscape_mutations.png
│       ├── permutation_test.png
│       ├── resistance_conservation_bars.png
│       └── conservation_by_group.png
```
