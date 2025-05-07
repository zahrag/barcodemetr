
import os
import pickle
from tabulate import tabulate
import numpy as np
from tqdm import tqdm
from pathlib import Path
from collections import defaultdict
import pandas as pd
from barcode_pwd_pandas import BarcodePWD as bar_pwd_panda
from barcode_pwd_spark import BarcodePWD as bar_pwd_spark


class BarcodeMetric:

    def __init__(self, metadata_file="", method="pandas", read_metadata=False):
        self.metadata = metadata_file
        self.df = self._read_metadata(metadata_file, read_metadata=False)
        self.taxonomy_ranks = ["phylum", "class", "order", "family", "subfamily", "genus", "species"]

        if method == "pandas":
            self.pwd = bar_pwd_panda()
        else:
            self.pwd = bar_pwd_spark()

        # Save all files in a pre-defined directory
        self.save_path = os.path.join(os.path.dirname(metadata_file), "barcode_analysis", method)
        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

    def _read_metadata(self, file, read_metadata=False):
        if not read_metadata:
            return None
        if file.endswith(".tsv"):
            return pd.read_csv(file, sep='\t', low_memory=False)
        elif file.endswith(".csv"):
            return pd.read_csv(file, low_memory=False)
        else:
            raise ValueError("Unsupported file extension. Use .tsv or .csv")

    def _sdi(self, sample_counts):
        proportions = np.array(sample_counts) / sum(sample_counts)
        sdi = -np.sum(proportions * np.log2(proportions + 1e-12))
        return sdi

    def build_hierarchy(self, df=None, taxonomy_ranks=None, path=None):
        """
        Build the barcode hierarchy, or load from a pickled file if path is provided.

        Args:
            df (pd.DataFrame, optional): DataFrame to use instead of reading from metadata.
            taxonomy_ranks (list, optional): Taxonomic ranks to process.
            path (str, optional): Path to a pickled hierarchy file to load.

        Returns:
            dict: Nested hierarchy of barcodes.
        """

        if path is not None and os.path.isfile(path):
            return self.open_pickle(path, enabled=True)

        df = df if df is not None else self.df
        taxonomy_ranks = taxonomy_ranks if taxonomy_ranks is not None else self.taxonomy_ranks

        data_hierarchy = defaultdict(lambda: defaultdict(lambda: {'barcodes': {}, 'sdi': 0.0}))
        for _, row in tqdm(df.iterrows(), total=len(df), desc="Building hierarchy"):
            barcode = row['dna_barcode']
            sample = row['processid']

            for rank in taxonomy_ranks:
                subgroup = row[rank]
                if pd.notna(subgroup):
                    subgroup_entry = data_hierarchy[rank][subgroup]
                    if barcode not in subgroup_entry['barcodes']:
                        subgroup_entry['barcodes'][barcode] = {'samples': []}
                    subgroup_entry['barcodes'][barcode]['samples'].append(sample)

        return data_hierarchy

    def compute_sdi(self, data_hierarchy, enabled=False):
        """
        Computes Shannon Diversity Index (SDI) for each (rank, subgroup) in the hierarchy.
        Updates the hierarchy in-place by storing 'sdi' per subgroup.
        SDI measures the diversity of barcodes within a subgroup.
        """

        if not enabled:
            return data_hierarchy

        for rank in data_hierarchy:
            for subgroup, subgroup_entry in tqdm(data_hierarchy[rank].items(), total=len(data_hierarchy[rank]), desc=f"Computing SDI of {rank}"):
                sample_counts = [
                    len(barcode_info['samples'])
                    for barcode_info in subgroup_entry['barcodes'].values()
                ]
                subgroup_entry['sdi'] = self._sdi(sample_counts)

        return self.convert_to_regular_dict(data_hierarchy)

    def compute_barcodes_statistics(self, data_hierarchy):

        barcode_stats = {}
        for rank in data_hierarchy:

            dna_unique_group_lst = []
            sdi_lst = []
            for subgroup, subgroup_entry in tqdm(data_hierarchy[rank].items(),
                                                 total=len(data_hierarchy[rank]),
                                                 desc=f"Barcode statistics of {rank}"):

                num_barcodes = len(subgroup_entry['barcodes'])
                sdi = subgroup_entry.get('sdi', 0.0)

                dna_unique_group_lst.append(num_barcodes)
                sdi_lst.append(sdi)

            barcode_stats[rank] = {
                'Average Distinct Barcodes': np.mean(dna_unique_group_lst),
                'Std.Dev Distinct Barcodes': np.std(dna_unique_group_lst),
                'Median Distinct Barcodes': np.median(dna_unique_group_lst),
                'Average Shannon Diversity Index': np.mean(sdi_lst),
                'Std.Dev Shannon Diversity Index': np.std(sdi_lst),
            }

        print(f"Identical DNA Barcode Statistics: {barcode_stats}")
        return barcode_stats

    def compute_pwd(self, data_hierarchy, taxonomy_ranks=None):
        """ Compute Damerau-Levenshtein pairwise distances of the identical DNA barcodes for a taxonomy rank"""
        taxonomy_ranks = taxonomy_ranks if taxonomy_ranks is not None else self.taxonomy_ranks
        for rank in taxonomy_ranks:
            self.pwd._rank_dist(data_hierarchy[rank], rank, path=self.save_path)

    def compute_pwd_statistics(self, data_hierarchy, taxonomy_ranks=None):

        taxonomy_ranks = taxonomy_ranks if taxonomy_ranks is not None else self.taxonomy_ranks

        pwd_stats = {}
        for rank in taxonomy_ranks:
            if rank in data_hierarchy:
                distances_root = f"{self.save_path}/distances/{rank}"
                if not os.path.exists(distances_root):
                    raise ValueError(f"The directory of distances files of {rank} does NOT exist.\n "
                                     f"Compute pairwise distances across {rank} first."
                                     )
                pwd_stats[rank] = self.pwd._rank_dist_stats(data_hierarchy[rank],
                                                            rank,
                                                            distances_root=distances_root)

        print(f"Identical DNA Barcode Pairwise Distance Statistics: {pwd_stats}")
        return pwd_stats

    def compute_full_statistics(self, data_hierarchy):

        barcode_stats = self.compute_barcodes_statistics(data_hierarchy)

        pwd_stats = self.compute_pwd_statistics(data_hierarchy)

        s_dict = defaultdict(list)
        ranks = list(data_hierarchy.keys())

        for rank in ranks:

            s_dict['Barcode Statistics'].append(rank)

            # Merge per-rank statistics
            merged_stats = {**barcode_stats.get(rank, {}), **pwd_stats.get(rank, {})}

            for key, value in merged_stats.items():
                s_dict[key].append(round(value, 4))

            print(f'Full Statistics of {rank} processed and collected.')

        return s_dict

    @staticmethod
    def convert_to_regular_dict(d):
        if isinstance(d, defaultdict):
            d = {k: BarcodeMetric.convert_to_regular_dict(v) for k, v in d.items()}
        elif isinstance(d, dict):
            d = {k: BarcodeMetric.convert_to_regular_dict(v) for k, v in d.items()}
        return d

    @staticmethod
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

    @staticmethod
    def save_in_pandas(df, panda_file, enabled=False):
        """Save distances as Pandas dataframe"""
        if not enabled:
            return
        df.reset_index(inplace=True, drop=True)
        if panda_file.endswith(".tsv"):
            df.to_csv(panda_file, sep='\t', index=False)
        elif panda_file.endswith(".csv"):
            df.to_csv(panda_file, index=False)
        else:
            raise ValueError("Unsupported file extension. Use .tsv or .csv")

    @staticmethod
    def create_pickle(data, pickle_file, enabled=False):
        if not enabled:
            return
        with open(pickle_file, 'wb') as f:
            pickle.dump(data, f)

    @staticmethod
    def open_pickle(pickle_file, enabled=False):
        if not enabled:
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






