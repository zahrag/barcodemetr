
import os

from tqdm import tqdm
import random
import textdistance

import pandas as pd
from barcode_alignment import perform_mafft_alignment
from pathlib import Path

# Get the path of the currently running script
current_directory = Path(__file__).parent

def read_tsv(tsv_file):
    df = pd.read_csv(tsv_file, sep='\t', low_memory=False)
    return df

class BarcodePWD(object):

    def _damerau_levenshtein_distance(self, aligned_sequences):
        """
        This function implements Damerau-Levenshtein distance.

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

    def _subgroup_distances(self, item, max_barcodes=1000):

        name, barcodes = item

        # Random sampling
        if (max_barcodes > 0) and (len(barcodes) > max_barcodes):
            sequences = random.sample(barcodes, max_barcodes)
        else:   # species
            sequences = barcodes

        # Perform alignment
        aligned_sequences = perform_mafft_alignment(sequences, name)

        # Compute pairwise distances
        distances = self._damerau_levenshtein_distance(aligned_sequences)

        return (name, distances)

    def _rank_dist(self, rank_hierarchy, rank, min_barcodes=4, max_barcodes=1000, chunk_size=1000):
        """
        This function process distance computation per taxonomic rank using pandas.
        :param rank_hierarchy: Data hierarchy at a specified taxonomic level.
        :param rank: Taxonomic group level (e.g., family, genus, species).
        :param min_barcodes: Minimum number of barcodes per rank to compute pairwise distances.
        :param max_barcodes: Maximum number of barcodes per rank randomly sampled to compute pairwise distances.
        :param chunk_size: Chunk size of the subgroups of the rank.
        """

        # Convert the dictionary to a list of tuples [(species, sequences), ...]
        print(f'Processing DNA barcodes pairwise distances across {rank} ...')

        tuple_list = [
            (subgroup, list(subgroup_entry['barcodes'].keys()))
            for subgroup, subgroup_entry in rank_hierarchy.items()
            if len(subgroup_entry['barcodes']) > min_barcodes
        ]
        # Create chunks of subgroups
        tuple_chunks = [tuple_list[i:i + chunk_size]
                        for i in range(0, len(tuple_list), chunk_size)
                        ]

        max_barcodes = 0 if rank == 'species' else max_barcodes     # Process all barcodes of the species
        for chk, chunk in enumerate(tuple_chunks):

            final_distances = None
            for cnt, item in tqdm(enumerate(chunk), total=len(chunk), desc="Processing groups"):

                group_name, df = self._subgroup_distances(item, max_barcodes=max_barcodes)

                df['group_name'] = group_name
                df['distance'] = df['distance'].astype(float)

                # Deselect DNA barcodes to save space
                final_distances = df[['distance', 'group_name']].copy() if final_distances is None \
                    else pd.concat([final_distances, df[['distance', 'group_name']]], ignore_index=True)

            # ---- Save dataframe of pairwise distance of subgroups in the chunk chk
            distances_dir = os.path.join(current_directory, 'distances', rank)
            if not os.path.exists(distances_dir):
                os.makedirs(distances_dir)
            final_distances.reset_index(inplace=True, drop=True)
            final_distances.to_csv(os.path.join(distances_dir, f'barcode_pwd_{rank}_chunk_{chk}.tsv'), sep='\t', index=False)

    def _rank_dist_stats(self, rank):

        """
        Compute pairwise distance statistics across taxonomic levels.
        :param rank: Taxonomic group level (e.g., family, genus, species).
        """

        distances_dir = os.path.join(current_directory, 'distances', rank)
        chks = self.extract_chunks(rank, distances_dir)

        print(f'Processing DNA barcodes statistics across {rank} ...')

        rank_stats = None
        df_pandas = None
        for chk_num, chk in tqdm(enumerate(chks), total=len(chks), desc="Processing statistics"):

            distances_tsv = os.path.join(distances_dir, f'barcodes_pwd_{rank}_chunk_{chk}.tsv')
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

    @staticmethod
    def extract_chunks(rank, dists_dir):
        dists_dir = Path(dists_dir)

        if not dists_dir.exists():
            raise ValueError(f"The directory of distances files of {rank} does NOT exist.\n"
                             f"Compute pairwise distances across {rank} first.")

        chk_values = [
            int(f.name.split("_chunk_")[1].split(".")[0])
            for f in dists_dir.iterdir()
            if f.name.startswith(f"barcode_pwd_{rank}_chunk_") and f.name.endswith(".tsv")
        ]
        return chk_values



