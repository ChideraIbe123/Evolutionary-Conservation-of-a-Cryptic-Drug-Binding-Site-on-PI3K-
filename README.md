# Evolutionary Conservation of a Cryptic Drug Binding Site on PI3Kα

**Chidera Ibe** — coibe2@illinois.edu

A bioinformatics analysis testing whether evolutionary conservation predicts which residues develop resistance mutations in patients treated with the allosteric PI3Kα inhibitor RLY-2608. All algorithms (Needleman-Wunsch with affine gaps, progressive multiple sequence alignment, neighbor-joining phylogeny, Shannon-entropy conservation scoring) are implemented from scratch in NumPy.

For more info, see [REPORT.md](REPORT.md).

## Project Structure

```
src/
  fetchers.py        Ensembl / UniProt / RCSB data fetching with caching
  utils.py           BLOSUM62 loader, FASTA I/O
  alignment.py       Needleman-Wunsch, profile alignment, progressive MSA, UPGMA
  tree.py            Neighbor-joining
  conservation.py    Shannon entropy conservation scoring

notebooks/
  01_alignment.ipynb               Pairwise alignment and progressive MSA
  02_phylogenetic_tree.ipynb       Neighbor-joining phylogeny
  03_conservation_structure.ipynb  Conservation scoring and PDB pocket mapping
  04_resistance_analysis.ipynb     Resistance mutation statistical analysis

data/
  raw/         Cached downloads from Ensembl, UniProt, RCSB
  processed/   MSAs (FASTA) and trees (Newick)
  results/     Conservation CSVs, summary tables, figures
```

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
jupyter lab notebooks/
```

## Run Order

1. `01_alignment.ipynb` — fetches sequences, runs pairwise NW alignments, builds progressive MSA
2. `02_phylogenetic_tree.ipynb` — builds neighbor-joining tree (also serves as guide tree for MSA)
3. `03_conservation_structure.ipynb` — computes conservation scores, identifies pocket residues from PDB 8TSD, maps onto 3D structure with py3Dmol
4. `04_resistance_analysis.ipynb` — statistical comparison of resistance mutation positions vs. conservation scores

## Key Findings

- PI3Kα is under extreme purifying selection: 87% of residues have conservation > 0.9
- The cryptic drug pocket is more conserved than the protein average (0.985 vs 0.969)
- All clinical resistance mutation positions (W780, Q859, E726, I817) are perfectly conserved across vertebrates
- Conservation does not discriminate resistance-prone sites because the entire protein is too conserved
- Resistance at hyper-conserved positions implies tumors pay a fitness cost to escape RLY-2608

## Data Sources

- PI3Kα orthologs: Ensembl REST API
- PI3K paralogs (α, β, γ, δ): UniProt (P42336, P42338, P48736, O00329)
- Co-crystal structure: PDB 8TSD
- Resistance mutations: Varkaris et al., *Cancer Discovery* 2024
