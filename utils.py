
from pathlib import Path
import pickle
from tabulate import tabulate
from collections import defaultdict


def convert_to_regular_dict(d):
    if isinstance(d, defaultdict):
        d = {k: convert_to_regular_dict(v) for k, v in d.items()}
    elif isinstance(d, dict):
        d = {k: convert_to_regular_dict(v) for k, v in d.items()}
    return d

def print_table(data_dict, title="", display_table=False):
    if not display_table:
        return
    print("\n\n" + "+" + "-" * 143 + "+")
    print(f"\t\t\t\t\t{title}")
    headings = list(data_dict.keys())
    values = list(data_dict.values())
    rows = zip(*values)
    formatted_rows = [
        [f'{val: .2f}' if isinstance(val, float) else str(val) for val in row]
        for row in rows
    ]
    print(tabulate(formatted_rows, headers=headings, tablefmt="grid"))

def save_in_pandas(df, path, _save=False):
    """Save data as Pandas dataframe"""
    df.reset_index(inplace=True, drop=True)
    if path.endswith(".tsv"):
        df.to_csv(path, sep='\t', index=False)
    elif path.endswith(".csv"):
        df.to_csv(path, index=False)
    else:
        raise ValueError("Unsupported file extension. Use .tsv or .csv")


def create_pickle(data=None, pickle_file=""):
    if not data:
        return
    with open(pickle_file, 'wb') as f:
        pickle.dump(data, f)

def open_pickle(pickle_file=""):
    if not pickle_file:
        return
    objects = []
    with (open(pickle_file, "rb")) as openfile:
        while True:
            try:
                objects.append(pickle.load(openfile))
            except EOFError:
                break
    results = objects[0]

    return results


def extract_chunks(rank, dists_dir):
    dists_dir = Path(dists_dir)

    if not dists_dir.exists():
        raise ValueError(f"The directory of distances files of {rank} does NOT exist.\n"
                         f"Compute pairwise distances across {rank} first.")

    chk_values = [
        int(f.name.split("_chunk_")[1].split(".")[0])
        for f in dists_dir.iterdir()
        if f.name.startswith(f"barcode_pwd_{rank}_chunk_") and f.name.endswith(".csv")
    ]
    return chk_values




