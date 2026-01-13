import gemmi

from .monlib import MonLib
from .reflections import DataItem
from .structure import remove_isolated_fragments
from .validation import validate


def prune(
    structure: gemmi.Structure,
    fphi_best: DataItem,
    fphi_diff: DataItem,
    fphi_calc: DataItem,
    residues: bool = True,
    chain_threshold: float = -2,
    residue_threshold: float = -5,
    monlib: MonLib = None,
) -> gemmi.Structure:
    print("Performing validation for pruning", flush=True)
    structure = structure.clone()
    monlib = monlib or MonLib(structure[0].get_all_residue_names())
    metrics = validate(structure, fphi_best, fphi_diff, fphi_calc, monlib)

    max_deleted = int(len(metrics) * 0.2)
    num_deleted = 0
    grouped = metrics.groupby("Chain")
    means = grouped.mean(numeric_only=True)
    for chain_name in means.sort_values("Score").index:
        score = means.loc[chain_name, "Score"]
        count = grouped.size().loc[chain_name]
        print(
            f"Chain {chain_name} has a score of {score} over {count} residues",
            flush=True,
        )
        if (
            means.loc[chain_name, "Score"] < chain_threshold
            and count <= 20
            and num_deleted + count <= max_deleted
        ):
            print("Deleting chain", chain_name, flush=True)
            del structure[0][chain_name]
            num_deleted += num_deleted
            metrics = metrics[metrics["Chain"] != chain_name]

    if not residues:
        return structure if num_deleted > 0 else None

    max_deleted = int(len(metrics) * 0.2)
    metrics = metrics[metrics["Score"] < residue_threshold]
    metrics.sort_values("Score", inplace=True)
    metrics = metrics.head(max_deleted)
    if len(metrics) == 0:
        return structure if num_deleted > 0 else None

    print(
        f"Deleting {len(metrics)} residues with scores < {residue_threshold}",
        flush=True,
    )
    to_delete = {(row["Chain"], row["SeqId"]) for _, row in metrics.iterrows()}
    for chain in structure[0]:
        for i, residue in reversed(list(enumerate(chain))):
            if (chain.name, str(residue.seqid)) in to_delete:
                del chain[i]

    print("Removing isolated residues (if any)", flush=True)
    for chain in structure[0]:
        remove_isolated_fragments(chain, monlib, max_length=1)

    return structure
