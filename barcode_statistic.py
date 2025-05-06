from tabulate import tabulate
import numpy as np
from tqdm import tqdm
from pathlib import Path
from collections import defaultdict
import pandas as pd
from barcode_pwd_pandas import BarcodePWD as bar_PWD_panda
from barcode_pwd_spark import BarcodePWD as bar_PWD_spark


# Get the path of the currently running script
current_directory = Path(__file__).parent

def convert_to_regular_dict(d):
    if isinstance(d, defaultdict):
        d = {k: convert_to_regular_dict(v) for k, v in d.items()}
    elif isinstance(d, dict):
        d = {k: convert_to_regular_dict(v) for k, v in d.items()}
    return d

class BarcodeMetric:

    def __init__(self, metadata_file="", method="pandas"):
        self.metadata = metadata_file
        self.df = self._read_metadata()
        self.taxonomy_ranks = ["dna_bin", "species", "genus", "subfamily", "family", "order", "class", "phylum"]

        if method == "pandas":
            self.pwd = bar_PWD_panda()
        else:
            self.pwd = bar_PWD_spark()


    def _read_metadata(self):
        return pd.read_csv(self.metadata, low_memory=False)

    def _generate_hierarchy(self, df=None, taxonomy_ranks=None):

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

    def _convert(self, default_dict):
        return convert_to_regular_dict(default_dict)

    def _sdi(self, sample_counts):
        proportions = np.array(sample_counts) / sum(sample_counts)
        sdi = -np.sum(proportions * np.log2(proportions + 1e-12))
        return sdi

    def compute_sdi(self, data_hierarchy):
        """
        Computes Shannon Diversity Index (SDI) for each (rank, subgroup) in the hierarchy.
        Updates the hierarchy in-place by storing 'sdi' per subgroup.
        SDI measures the diversity of barcodes within a subgroup.
        """
        for rank in data_hierarchy:
            for subgroup, subgroup_entry in data_hierarchy[rank].items():
                sample_counts = [
                    len(barcode_info['samples'])
                    for barcode_info in subgroup_entry['barcodes'].values()
                ]
                subgroup_entry['sdi'] = self._sdi(sample_counts)
        return self._convert(data_hierarchy)

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
        return barcode_stats

    def compute_pwd(self, data_hierarchy, taxonomy_ranks=None):
        """ Compute Damerau-Levenshtein pairwise distances of the identical DNA barcodes for a taxonomy rank"""
        taxonomy_ranks = taxonomy_ranks if taxonomy_ranks is not None else self.taxonomy_ranks
        for rank in taxonomy_ranks:
            self.pwd._rank_dist(data_hierarchy, rank)

    def compute_pwd_statistics(self, data_hierarchy, ranks=None):
        if ranks is None:
            ranks = list(data_hierarchy.keys())

        pwd_stats = {}
        for rank in ranks:
            if rank in data_hierarchy:
                pwd_stats[rank] = self.pwd._rank_dist_stats(data_hierarchy[rank], rank)
        return pwd_stats

    def compute_full_statistics(self, data_hierarchy):

        barcode_stats = self.compute_barcodes_statistics(data_hierarchy)

        pwd_stats = self.compute_pwd_statistics(data_hierarchy)

        s_dict = defaultdict(list)
        ranks = list(data_hierarchy.keys())

        for rank in ranks:
            s_dict['DNA Statistics'].append(rank)

            # Merge per-rank statistics
            merged_stats = {**barcode_stats.get(rank, {}), **pwd_stats.get(rank, {})}

            for key, value in merged_stats.items():
                s_dict[key].append(round(value, 4))

            print(f'Statistics of {rank} processed and collected.')
        self.print_table(s_dict, "Barcode Statistics", print_table=True)

        return s_dict

    def print_table(self, data_dict, title, print_table=False):

        if not print_table:
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


def main():
    barmetr = BarcodeMetric(metadata="", method="spark")
    hierarchy = barmetr._generate_hierarchy()
    hierarchy = barmetr.compute_sdi(hierarchy)

    barmetr.compute_pwd(hierarchy, rank="species")


if __name__ == "__main__":
    main()





