
import os
import numpy as np
from tqdm import tqdm

from utils import *
from barcode_pwd_pandas import BarcodePWD as bar_pwd_panda
from barcode_pwd_spark import BarcodePWD as bar_pwd_spark



class BarcodeMetric:
    """ Main class to compute barcode statistics."""

    def __init__(self, method="pandas", max_barcode_length=625, metadata_file="", save_path="", load_metadata=False):
        """
        :param method (str): Method to compute barcode statistics; spark or pandas.
        :param max_barcode_length (int, Optional): Maximum barcode length applied in alignment.
        :param metadata_file (str): Path to metadata file.
        :param save_path (str): Path to save resulting files.
        :param load_metadata (boolean): If load the metadata file; True only to compute pairwise distances.
        """

        self.metadata = metadata_file
        self.df = load_from_pandas(metadata_file, load_file=load_metadata)
        self.taxonomy_ranks = ["phylum", "class", "order", "family", "subfamily", "genus", "species"]

        # Save all files in a pre-defined directory
        self.save_path = save_path
        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

        if method == "pandas":
            self.pwd = bar_pwd_panda(save_path=self.save_path, max_barcode_length=max_barcode_length)
        else:
            self.pwd = bar_pwd_spark(save_path=self.save_path, max_barcode_length=max_barcode_length)

    def sdi(self, sample_counts):
        """ Compute the Shannon Diversity Index across a subgroup from its counts of unique DNA barcodes"""
        proportions = np.array(sample_counts) / sum(sample_counts)
        sdi = -np.sum(proportions * np.log2(proportions + 1e-12))
        return sdi

    def build_hierarchy(self, df=None, taxonomy_ranks=None, path=None):
        """
        Build the barcode hierarchy, or load from a pickled file if path is provided.

        Args:
            df (pd.DataFrame, optional): DataFrame to use instead of reading from file.
            taxonomy_ranks (list, optional): Taxonomic ranks to process.
            path (str): Path to a pickled hierarchy file to load.

        Returns:
            dict: Nested hierarchy of barcodes.
        """

        if path is not None and os.path.isfile(path):
            return open_pickle(pickle_file=path)

        df = df if df is not None else self.df
        taxonomy_ranks = taxonomy_ranks if taxonomy_ranks is not None else self.taxonomy_ranks

        ranked_data = defaultdict(lambda: defaultdict(lambda: {'barcodes': {}, 'sdi': 0.0}))
        for _, row in tqdm(df.iterrows(), total=len(df), desc="Building hierarchy"):
            barcode = row['dna_barcode']
            sample = row['processid']

            for rank in taxonomy_ranks:
                subgroup = row[rank]
                if pd.notna(subgroup):
                    subgroup_entry = ranked_data[rank][subgroup]
                    if barcode not in subgroup_entry['barcodes']:
                        subgroup_entry['barcodes'][barcode] = {'samples': []}
                    subgroup_entry['barcodes'][barcode]['samples'].append(sample)

        return ranked_data

    def compute_sdi(self, ranked_data, path=None):
        """
        Computes Shannon Diversity Index (SDI) across subgroups of each rank in the hierarchy.
        Updates the ranked data dict in-place by storing 'sdi' per subgroup.
        SDI measures the diversity of barcodes within a subgroup.
        """

        if path is not None and os.path.isfile(path):
            return ranked_data

        for rank in ranked_data:
            for subgroup, subgroup_entry in tqdm(ranked_data[rank].items(), total=len(ranked_data[rank]), desc=f"Computing SDI of {rank}"):
                sample_counts = [
                    len(barcode_info['samples'])
                    for barcode_info in subgroup_entry['barcodes'].values()
                ]
                subgroup_entry['sdi'] = self.sdi(sample_counts)

        return convert_to_regular_dict(ranked_data)

    def compute_barcodes_statistics(self, ranked_data):
        """Compute unique DNA barcode descriptive statistics"""

        barcode_stats = {}
        for rank in ranked_data:

            dna_unique_group_lst = []
            sdi_lst = []
            for subgroup, subgroup_entry in tqdm(ranked_data[rank].items(),
                                                 total=len(ranked_data[rank]),
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

        print(f"\nUnique DNA Barcode Statistics: {barcode_stats}\n")
        return barcode_stats

    def compute_pwd(self, ranked_data, taxonomy_ranks=None, _enabled=False):
        """ Compute Damerau-Levenshtein pairwise distances of the unique DNA barcodes for a taxonomy rank"""
        if not _enabled:
            return
        taxonomy_ranks = taxonomy_ranks if taxonomy_ranks is not None else self.taxonomy_ranks
        for rank in taxonomy_ranks:
            self.pwd._rank_dists(ranked_data[rank], rank, path=self.save_path)
        self.pwd._end_spark()

    def compute_pwd_statistics(self, ranked_data, taxonomy_ranks=None, save_distances_pandas=False):
        """ Compute Damerau-Levenshtein pairwise distances statistics"""
        taxonomy_ranks = taxonomy_ranks if taxonomy_ranks is not None else self.taxonomy_ranks

        pwd_stats = {}
        for rank in taxonomy_ranks:
            if rank in ranked_data:

                distances_root = os.path.join(self.save_path, rank)
                if not os.path.exists(distances_root):
                    raise FileNotFoundError(
                        f"Distance files for rank '{rank}' not found at: {distances_root}.\n"
                        "Please compute pairwise distances before proceeding."
                    )

                chunks = extract_chunks(rank, distances_root, method="spark")
                if not chunks:
                    raise ValueError(
                        f"No distance chunks found for rank '{rank}' in: {distances_root}.\n"
                        "Ensure pairwise distances were computed and saved correctly."
                    )

                print(f"\nFound {len(chunks)} saved chunks for rank '{rank}':\n{chunks}\n")

                pwd_stats[rank] = self.pwd._rank_dist_stats(rank, chunks,
                                                            distances_root=distances_root,
                                                            save_distances_pandas=save_distances_pandas,
                                                            )

        self.pwd._end_spark()
        print(f"Unique DNA Barcode Pairwise Distance Statistics: {pwd_stats}")
        return pwd_stats

    def compute_full_statistics(self, ranked_data, save_distances_pandas=False, _enabled=False):
        """ Compute full statistics of the unique DNA barcodes"""
        if not _enabled:
            return

        barcode_stats = self.compute_barcodes_statistics(ranked_data)

        pwd_stats = self.compute_pwd_statistics(ranked_data, save_distances_pandas=save_distances_pandas)

        s_dict = defaultdict(list)
        ranks = list(ranked_data.keys())

        for rank in ranks:

            s_dict['Barcode Statistics'].append(rank)

            # Merge per-rank statistics
            merged_stats = {**barcode_stats.get(rank, {}), **pwd_stats.get(rank, {})}

            for key, value in merged_stats.items():
                s_dict[key].append(round(value, 4))

            print(f'Full Statistics of {rank} Processed and Collected.')

        return s_dict

