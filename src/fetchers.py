"""
Data fetching utilities with local file caching.
Supports Ensembl REST API, UniProt, and RCSB PDB.
"""

import os
import time
import json
import requests

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")


def cached_fetch(url, cache_path, headers=None, binary=False):
    """Download a URL and cache the result locally. Returns file contents."""
    if os.path.exists(cache_path):
        mode = "rb" if binary else "r"
        with open(cache_path, mode) as f:
            return f.read()

    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    resp = requests.get(url, headers=headers, timeout=60)
    resp.raise_for_status()

    mode = "wb" if binary else "w"
    content = resp.content if binary else resp.text
    with open(cache_path, mode) as f:
        f.write(content)
    return content


def fetch_ensembl_orthologs(gene_symbol="PIK3CA", species="human"):
    """
    Fetch ortholog list for a gene from Ensembl REST API.
    Returns list of dicts with keys: species, protein_id, perc_id, type.
    """
    cache_path = os.path.join(RAW_DIR, "orthologs", f"{gene_symbol}_orthologs.json")
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            return json.load(f)

    url = (
        f"https://rest.ensembl.org/homology/symbol/{species}/{gene_symbol}"
        f"?type=orthologues&format=condensed&content-type=application/json"
    )
    resp = requests.get(url, headers={"Content-Type": "application/json"}, timeout=120)
    resp.raise_for_status()
    data = resp.json()

    orthologs = []
    for homology in data["data"][0]["homologies"]:
        # Condensed format puts fields directly in the homology dict
        orthologs.append({
            "species": homology["species"],
            "protein_id": homology["protein_id"],
            "perc_id": homology.get("perc_id", 0),
            "type": homology.get("type", "ortholog_one2one"),
        })

    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(orthologs, f, indent=2)

    return orthologs


def fetch_ensembl_sequences(protein_ids, batch_size=50):
    """
    Fetch protein sequences from Ensembl in batches.
    Returns dict mapping protein_id -> sequence string.
    """
    cache_dir = os.path.join(RAW_DIR, "orthologs", "sequences")
    os.makedirs(cache_dir, exist_ok=True)

    sequences = {}

    # Load already-cached sequences
    uncached_ids = []
    for pid in protein_ids:
        cache_path = os.path.join(cache_dir, f"{pid}.fasta")
        if os.path.exists(cache_path):
            with open(cache_path) as f:
                content = f.read()
            # Parse simple FASTA
            lines = content.strip().split("\n")
            seq = "".join(l for l in lines if not l.startswith(">"))
            sequences[pid] = seq
        else:
            uncached_ids.append(pid)

    # Batch-fetch uncached sequences
    for i in range(0, len(uncached_ids), batch_size):
        batch = uncached_ids[i:i + batch_size]
        url = "https://rest.ensembl.org/sequence/id"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        body = {"ids": batch}

        resp = requests.post(url, json=body, headers=headers, timeout=120)
        resp.raise_for_status()
        results = resp.json()

        for entry in results:
            pid = entry["id"]
            seq = entry["seq"]
            sequences[pid] = seq

            # Cache individually
            cache_path = os.path.join(cache_dir, f"{pid}.fasta")
            with open(cache_path, "w") as f:
                f.write(f">{pid}\n{seq}\n")

        # Rate limiting
        if i + batch_size < len(uncached_ids):
            time.sleep(1)

    return sequences


# UniProt accessions for human PI3K isoforms
PI3K_PARALOGS = {
    "PI3Ka_PIK3CA": "P42336",
    "PI3Kb_PIK3CB": "P42338",
    "PI3Kg_PIK3CG": "P48736",
    "PI3Kd_PIK3CD": "O00329",
}


def fetch_uniprot_sequence(accession, name=None):
    """Fetch a protein sequence from UniProt in FASTA format."""
    if name is None:
        name = accession
    cache_path = os.path.join(RAW_DIR, "paralogs", f"{name}.fasta")
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    content = cached_fetch(url, cache_path)

    # Parse FASTA
    lines = content.strip().split("\n")
    header = lines[0]
    seq = "".join(l for l in lines[1:] if not l.startswith(">"))
    return header, seq


def fetch_all_paralogs():
    """Fetch all 4 human PI3K paralog sequences. Returns dict name -> sequence."""
    paralogs = {}
    for name, accession in PI3K_PARALOGS.items():
        _, seq = fetch_uniprot_sequence(accession, name)
        paralogs[name] = seq
    return paralogs


def fetch_pdb(pdb_id="8TSD"):
    """Fetch a PDB file from RCSB. Returns the file contents as a string."""
    cache_path = os.path.join(RAW_DIR, "structures", f"{pdb_id}.pdb")
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    return cached_fetch(url, cache_path)


def fetch_human_pik3ca_sequence():
    """Convenience: fetch human PI3K-alpha sequence."""
    _, seq = fetch_uniprot_sequence("P42336", "PI3Ka_PIK3CA")
    return seq
