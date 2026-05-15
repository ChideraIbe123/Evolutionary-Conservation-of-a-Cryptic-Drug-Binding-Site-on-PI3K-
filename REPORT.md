# Evolutionary Conservation of a Cryptic Drug Binding Site on PI3Kα

**Chidera Ibe (coibe2)**

Code: <https://github.com/ChideraIbe123/Evolutionary-Conservation-of-a-Cryptic-Drug-Binding-Site-on-PI3K->

## Intro

PI3Kα is a lipid kinase that's mutated in about 40% of HR+ breast cancers, so it's a big drug target. The common driver mutations are H1047R, E545K, and E542K. Older drugs like alpelisib work but have nasty side effects (hyperglycemia, rash) because they also hit wild-type PI3Kα in healthy tissue. RLY-2608 (D.E. Shaw, 2024) is newer and more selective. It binds a cryptic allosteric pocket at the C2/kinase interface that mostly only forms when PI3Kα is mutated. Varkaris et al. (2024, *Cancer Discovery*) reported resistance mutations from patients on PI3Kα inhibitors at W780, Q859, E726, and I817, all near the pocket.

I wanted to check whether these resistance positions are evolutionarily different from the rest of the protein. The idea: if the pocket is super conserved, resistance should be hard (mutating a conserved residue should hurt the protein). If the pocket is variable, then maybe conservation could tell us which sites are likely to mutate. I tested this using Needleman-Wunsch with affine gaps, progressive MSA, and neighbor-joining, all implemented from scratch.

## Methods

### Data

I used PDB 8TSD (the co-crystal of PI3Kα + RLY-2608, ligand is XUZ). For orthologs I pulled 20 vertebrate PI3Kα sequences from Ensembl (primates, rodents, carnivores, ungulates, marsupials, birds, reptiles, amphibians, fish). For paralogs I grabbed the four human PI3K isoforms from UniProt (PIK3CA, PIK3CB, PIK3CG, PIK3CD). Resistance mutations are from the Varkaris paper.

### Needleman-Wunsch (affine gaps)

Wrote NW from scratch in NumPy with three DP matrices: $M$ for match/mismatch, $I_x$ for gap-in-y, and $I_y$ for gap-in-x. Used BLOSUM62 with gap-open $-10$ and gap-extend $-0.5$. The traceback was the annoying part. You have to handle the case where one index hits 0 before the other, which isn't obvious from the standard pseudocode. I checked my scores against Biopython's `PairwiseAligner` and they matched.

### Neighbor-joining

Implemented Saitou-Nei from scratch: get net divergence, build the Q-matrix, join the min-Q pair, get branch lengths, update distances, repeat. The pairwise distances came from NW % identities with the Jukes-Cantor correction. To sanity-check, I made sure the tree had the clades I expected: human/chimp together, mouse/rat together, chicken/turkey together, fish as an outgroup. All good.

### Progressive MSA

Used the NJ tree as the guide tree and did profile-profile alignment at each internal node. Each merge is basically NW where the column score between two profile columns is the frequency-weighted average of BLOSUM62 scores over all amino acid pairs. Ortholog MSA ended up at 20 sequences × 1073 columns. Paralog MSA is 4 sequences × 1206 columns.

### Conservation

Per-column Shannon entropy (ignoring gaps), normalized so 1 = fully conserved and 0 = random:

$$
H = -\sum_i p_i \log_2 p_i, \qquad C = 1 - \frac{H}{\log_2 20}.
$$

Then I mapped column indices back to human PI3Kα residue numbers using the non-gap positions of the human row.

### Finding the pocket

Parsed 8TSD with Biopython (chain A only). Used `NeighborSearch` to grab every residue with at least one atom within 5 Å of the ligand. Stuck the conservation scores in the B-factor column and viewed it in py3Dmol (cartoon colored by conservation, pocket as sticks, ligand in yellow).

### Stats

Mann-Whitney U on resistance vs non-resistance conservation scores, a permutation test (10k iterations) to compare resistance positions to random residue sets, and Fisher's exact on resistance × pocket. All using SciPy.

## Results

### PI3Kα barely changes across vertebrates

87.3% of the 1068 residues have conservation $> 0.9$. Median is 1.0, mean is 0.969. Basically the whole protein is conserved.

### The pocket is more conserved than the rest

`NeighborSearch` found 26 residues within 5 Å of RLY-2608 (positions 808 to 1025). Mean conservation in the pocket is 0.985 vs 0.969 elsewhere. Makes sense, since a lot of these residues are at the C2/kinase interface which probably matters even when the pocket is closed.

### Every resistance site is 100% conserved

W780, Q859, E726, I817 all score 1.0. So do the primary driver positions H1047, E545, E542. Every cancer-related mutation site in this list is invariant across all 20 species.

### Conservation can't really predict resistance though

Since the resistance positions *and* most of the protein are pinned at 1.0, Mann-Whitney U comes back at $p = 0.28$, permutation test at $p = 0.36$. Not significant. The problem isn't really that conservation is wrong, it's that there's no signal to find when almost every position is conserved. To actually predict resistance you'd probably need to add structural features (ligand contact, solvent accessibility, etc.).

### Paralogs explain selectivity

The four human PI3K paralogs only share 29 to 56% identity pairwise (α/β 38.6%, α/δ 38.9%, β/δ 56.4%, γ around 29% from everything). Out of the 26 pocket residues in PI3Kα, only 15 are the same across all four. The other 11 differ in at least one isoform. That's probably why RLY-2608 is selective for α, since the pocket looks different in β/γ/δ.

## Discussion

The interesting thing is that resistance mutations happen at the *most* conserved positions, not the variable ones. I was kind of expecting the opposite, that resistance would show up at positions evolution had let drift, since those should be the easiest to mutate without breaking the protein. Instead, the drug is pushing tumors to mutate residues that don't vary across the 20 species I sampled.

I think this is actually good news for the drug. If a tumor has to mutate something like W780 to escape RLY-2608, it's probably paying a real cost (lost activity, broken regulation, less stability). The pocket being super conserved also means there aren't that many positions that could even tolerate a resistance mutation, so the tumor's options are limited. And the fact that the pocket is conserved across α orthologs but different across paralogs is why the drug is selective in the first place.

Limitations: only 20 vertebrate species, and for a protein this conserved that might not be enough variation to see anything subtle. Adding invertebrate or fungal homologs would help but would also be harder to align. The resistance list is also pretty small (5 mutations, 4 positions) so the stats are weak. And conservation is just one signal. Combining it with structural stuff in some kind of model would probably do better.

## Conclusion

I implemented NW with affine gaps, progressive MSA, and neighbor-joining from scratch, applied them to 20 vertebrate PI3Kα orthologs and the 4 human paralogs, and mapped conservation onto the RLY-2608 binding site. The pocket is more conserved than the protein average, every clinical resistance site is 100% conserved, and paralog divergence at the pocket explains the drug's α-selectivity. Conservation by itself can't predict resistance because the protein is just too conserved, but the fact that resistance happens at fully-conserved sites at all suggests tumors are paying a fitness cost, which is good for the drug.
