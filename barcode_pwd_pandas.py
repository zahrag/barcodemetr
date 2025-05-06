
import os

from tqdm import tqdm
import random
import textdistance

import pandas as pd
from barcode_alignment import perform_mafft_alignment
from pathlib import Path

# Get the path of the currently running script
current_directory = Path(__file__).parent

def make_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)

def make_tsv(df, name=None, path=None):
    df_ = pd.DataFrame(df)
    df_.reset_index(inplace=True, drop=True)
    df_.to_csv(os.path.join(path, name), sep='\t', index=False)

def read_tsv(tsv_file):
    df = pd.read_csv(tsv_file, sep='\t', low_memory=False)
    return df

class BarcodePWD(object):
    """ Pairwise distance analysis of identical DNA barcode sequences using 'Damerau-Levenshtein Distance' """
    def __init__(self):

        self.min_seq_num = 3
        self.group_chunk_size = 1000
        self.seq_chunk_size = 1000

    def _pwd(self, aligned_sequences):
        """
        This function implements Damerau-Levenshtein distance using Pandas.

        Damerau-Levenshtein distance represents how many edits are needed to transform one sequence into another,
        where an edit can be an insertion, deletion, substitution, or transposition of adjacent characters (nucleotides).

        :param aligned_sequences: List of aligned sequences.
        :return: DataFrame with pairwise Levenshtein distances.
        """

        # Create a DataFrame with all pairwise combinations
        seq_df = pd.DataFrame(aligned_sequences, columns=["sequence"])
        cross_df = seq_df.merge(seq_df, how="cross", suffixes=("", "2"))

        # Calculate Damerau-Levenshtein distance
        cross_df["distance"] = cross_df.apply(
            lambda row: textdistance.damerau_levenshtein(row["sequence"], row["sequence2"]), axis=1
        )

        # Keep only unique pairs where sequence < sequence2
        unique_distances_df = cross_df[cross_df["sequence"] < cross_df["sequence2"]]

        return unique_distances_df[["sequence", "sequence2", "distance"]]

    def _subgroup_distances(self, item, max_len=1000):

        name, sequences = item

        # Random sampling
        if max_len > 0:
            sequences = random.sample(sequences, max_len) if len(sequences) > max_len else sequences

        # Perform alignment
        aligned_sequences = perform_mafft_alignment(sequences, name)

        # Compute pairwise distances
        distances = self._pwd(aligned_sequences)

        return (name, distances)

    def _rank_distances(self, data_dict, rank, process=False):
        """
        This function process distance computation per taxonomic rank using pandas.
        """

        if not process:
            return

        # Convert the dictionary to a list of tuples [(species, sequences), ...]
        print(f'Processing DNA barcodes pairwise distances across {rank} ...')

        # Process all DNA barcode sequences of the taxonomic rank (e.g., species)
        self.seq_chunk_size = 0 if rank == 'species' else self.seq_chunk_size

        tuple_list = [(group_name, data['unique_barcodes'])
                      for group_name, data in data_dict.items()
                      if len(data['unique_barcodes']) > self.min_seq_num ]

        tuple_chunks = [tuple_list[i:i + self.group_chunk_size] for i in range(0, len(tuple_list), self.group_chunk_size)]


        chk_ended = 0
        for chk, chunk in enumerate(tuple_chunks):

            if chk < chk_ended:
                continue

            final_distances = None
            for cnt, item in tqdm(enumerate(chunk), total=len(chunk), desc="Processing groups"):

                group_name, df = self._subgroup_distances(item, max_len=self.seq_chunk_size)

                df['group_name'] = group_name
                df['distance'] = df['distance'].astype(float)

                # Deselect DNA barcodes to save space
                final_distances = df[['distance', 'group_name']].copy() if final_distances is None \
                    else pd.concat([final_distances, df[['distance', 'group_name']]], ignore_index=True)

            # ---- Save the DataFrame to TSV format
            distances_dir = os.path.join(current_directory, 'distances', rank)
            make_directory(distances_dir)
            make_tsv(final_distances, name=f'barcode_pwd_{rank}_chunk_{chk}.tsv', path=distances_dir)


    def _rank_statistics(self, rank, max_chk=10, process=False):

        if not process:
            return None

        print(f'Processing DNA barcodes statistics across {rank} ...')

        rank_stats = None
        df_pandas = None
        for chk in tqdm(range(max_chk), total=max_chk, desc="Processing statistics"):

            distances_dir = os.path.join(current_directory, 'distances', rank)
            distances_tsv = os.path.join(distances_dir, f'barcode_pwd_{rank}_chunk_{chk}.tsv')
            if not os.path.exists(distances_dir):
                continue

            df = read_tsv(distances_tsv)
            df_pandas = df if df_pandas is None else pd.concat([df_pandas, df], ignore_index=True)

            # Group by 'group_name' and calculate statistics for each group of the chunk
            stats_df = df.groupby("group_name").agg(
                mean=("distance", "mean"),
                variance=("distance", "var"),
                stddev=("distance", "std"),
                min=("distance", "min"),
                max=("distance", "max")
            ).reset_index()

            # Combine the stats DataFrames
            rank_stats = stats_df if rank_stats is None else pd.concat([rank_stats, stats_df], ignore_index=True)

        aggregated_rank_stat = rank_stats.agg({
            "mean": "mean",
            "variance": "mean",
            "stddev": "mean",
            "min": "mean",
            "max": "mean"
        }).rename({
            "mean": "mean_of_means",
            "variance": "mean_of_variance",
            "stddev": "mean_of_stddev",
            "min": "mean_of_min",
            "max": "mean_of_max"
        })

        # Show the result
        print(aggregated_rank_stat)

        # Collect the results into a row
        aggregated_row = aggregated_rank_stat.to_frame().T

        # Convert the row to a dictionary
        rank_stats_dict = {
            "Mean": aggregated_row.loc[0, "mean_of_means"],
            "Variance": aggregated_row.loc[0, "mean_of_variance"],
            "Std.Dev": aggregated_row.loc[0, "mean_of_stddev"],
            "Min": aggregated_row.loc[0, "mean_of_min"],
            "Max": aggregated_row.loc[0, "mean_of_max"]
        }

        return rank_stats_dict


