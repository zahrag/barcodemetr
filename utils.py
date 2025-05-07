
from pathlib import Path
import pickle
from tabulate import tabulate
from collections import defaultdict
import pandas as pd


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

def save_in_pandas(data, path, _save=False):
    """Save data as Pandas dataframe"""
    if not _save:
        return

    # Convert based on data type
    if isinstance(data, dict) and all(isinstance(v, list) for v in data.values()):
        df = pd.DataFrame(data)
    elif hasattr(data, "toPandas"):  # likely a PySpark DataFrame
        df = data.toPandas()
    elif isinstance(data, pd.DataFrame):
        df = data
    else:
        raise TypeError("Unsupported data type for saving as pandas DataFrame")

    df.reset_index(inplace=True, drop=True)

    if path.endswith(".tsv"):
        df.to_csv(path, sep='\t', index=False)
    elif path.endswith(".csv"):
        df.to_csv(path, index=False)
    else:
        raise ValueError("Unsupported file extension. Use .tsv or .csv")


def load_from_pandas(file, load_file=False):
    if not load_file:
        return None

    if file.endswith(".tsv"):
        return pd.read_csv(file, sep='\t', low_memory=False)
    elif file.endswith(".csv"):
        return pd.read_csv(file, low_memory=False)
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


def extract_chunks(rank, dists_dir, method="spark"):

    dists_dir = Path(dists_dir)

    if not dists_dir.exists():
        raise ValueError(f"The directory of distances files of {rank} does NOT exist.\n"
                         f"Compute pairwise distances across {rank} first.")

    if method == "panda":

        chk_values = [
            int(f.name.split("_chunk_")[1].split(".")[0])
            for f in dists_dir.iterdir()
            if f.name.startswith(f"barcode_pwd_{rank}_chunk_") and f.name.endswith(".csv")
        ]

        return chk_values

    elif method == "spark":

        chk_values = [
            int(f.name.split("_")[1])
            for f in dists_dir.iterdir()
            if f.is_dir() and f.name.startswith("chunk_") and f.name.split("_")[1].isdigit()
        ]
        return chk_values

    else:
        raise ValueError(f"Unsupported method {method}")





