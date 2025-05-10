
from tqdm import tqdm
import random
import textdistance
from barcode_alignment import BarcodeAlignment

from utils import *

class BarcodePWD(object):
    """ DNA Barcodes Pairwise Distance Calculator by Pandas"""

    def __init__(self, max_barcode_length=None, save_path=None):
        """
        :param max_barcode_length: Maximum DNA barcode length applied for alignment.
        :param save_path: Path to save resulting files.
        """
        self.save_path = save_path
        # Initialize Barcode Alignment
        self.mafft_align = BarcodeAlignment(save_path=save_path, max_barcode_length=max_barcode_length)

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

    def _subgroup_dists(self, subgroup, max_barcodes=1000):
        """
        Compute pairwise distances of a subgroup.
        :param subgroup: The corresponding subgroup of a taxonomic rank.
        :param max_barcodes: Maximum number of barcodes randomly sampled in the subgroup.
        :return: Name and pairwise distances of the subgroup.
        """

        name, barcodes = subgroup

        # Random sampling
        if (max_barcodes > 0) and (len(barcodes) > max_barcodes):
            sequences = random.sample(barcodes, max_barcodes)
        else:   # species
            sequences = barcodes

        # Perform alignment
        aligned_sequences = self.mafft_align.perform_mafft_alignment(sequences, name)

        # Compute pairwise distances
        distances = self._damerau_levenshtein_distance(aligned_sequences)

        return (name, distances)

    def _rank_dist(self, rank_hierarchy, rank, min_barcodes=4, max_barcodes=1000, chunk_size=1000, path=None):
        """
        This function process distance computation per taxonomic rank using pandas.

        :param rank_hierarchy: Data hierarchy at a specified taxonomic level.
        :param rank: Taxonomic group level (e.g., family, genus, species).
        :param min_barcodes: Minimum number of barcodes per subgroups to compute pairwise distances.
        :param max_barcodes: Maximum number of barcodes per subgroups randomly sampled to compute pairwise distances.
        :param chunk_size: Chunk size of the subgroups of the rank.
        :param path: Path to save the resulting dataframes.
        """

        tuple_list = [
            (subgroup, list(subgroup_entry['barcodes'].keys()))
            for subgroup, subgroup_entry in rank_hierarchy.items()
            if len(subgroup_entry['barcodes']) > min_barcodes
        ]

        tuple_chunks = [tuple_list[i:i + chunk_size] for i in range(0, len(tuple_list), chunk_size) ]

        max_barcodes = 0 if rank == 'species' else max_barcodes     # Process all barcodes of the species
        for chk, chunk in enumerate(tuple_chunks):

            final_distances = None
            for cnt, item in tqdm(enumerate(chunk), total=len(chunk), desc=f"Processing groups of {rank}"):

                group_name, df = self._subgroup_dists(item, max_barcodes=max_barcodes)

                df['group_name'] = group_name
                df['distance'] = df['distance'].astype(float)

                # Deselect DNA barcodes to save space
                final_distances = df[['distance', 'group_name']].copy() if final_distances is None \
                    else pd.concat([final_distances, df[['distance', 'group_name']]], ignore_index=True)

            # ---- Save dataframe of pairwise distance of subgroups in the chunk chk
            path = path if path is not None else self.save_path

            file_path = os.path.join(path, rank, f"chunk_{chk}", f"barcodes_pwd_{rank}_chunk_{chk}.csv")
            if not os.path.exists(os.path.dirname(file_path)):
                os.makedirs(os.path.dirname(file_path))

            save_in_pandas(final_distances, file_path, _save=True)

    def _rank_dist_stats(self, rank, chunks, distances_root=None):

        """
        Compute pairwise distance statistics across taxonomic levels.

        :param rank: Taxonomic group level (e.g., family, genus, species).
        :param chunks: Available chunks of the subgroups of the rank.
        :param distances_root: Root path of the distances directory.
        """

        rank_stats = None
        df_pandas = None
        for chk_num, chk in tqdm(enumerate(chunks), total=len(chunks), desc=f"Processing statistics of {rank}"):

            file_path = os.path.join(distances_root, f"chunk_{chk}", f"barcodes_pwd_{rank}_chunk_{chk}.csv")
            df = pd.read_csv(file_path, low_memory=False)

            df_pandas = df if df_pandas is None else pd.concat([df_pandas, df], ignore_index=True)

            # Group by 'subgroup_name' and calculate statistics for each group of the chunk
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
        # print(aggregated_rank_stat)

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

