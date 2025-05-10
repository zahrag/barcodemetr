

import glob
from tqdm import tqdm
import random
from pyspark.sql import SparkSession
from pyspark.sql import functions as F
from pyspark.sql.functions import col
import textdistance

from barcode_alignment import BarcodeAlignment
from utils import *


class BarcodePWD(object):
    """ DNA Barcodes Pairwise Distance Calculator by Apache Spark"""

    def __init__(self, max_barcode_length=None, save_path=None):
        """
       :param max_barcode_length: Maximum DNA barcode length applied for alignment.
       :param save_path: Path to save resulting files.
       """
        self.save_path = save_path

        # try
        self.spark = self.initialize_spark()
        # Determine the number of available cores
        self.num_cores = self.spark.sparkContext.defaultParallelism

        # Set the number of partitions to 2-4 times the number of cores
        self.num_partitions = max(2, min(4 * self.num_cores, 100))  # Limit to a reasonable max, e.g., 100

        # Initialize Barcode Alignment
        self.mafft_align = BarcodeAlignment(save_path=save_path, max_barcode_length=max_barcode_length)

    def initialize_spark(self):
        return SparkSession.builder \
    .appName("Pairwise Levenshtein") \
    .config("spark.driver.memory", "10g") \
    .config("spark.executor.memory", "10g") \
    .config("spark.driver.memoryOverhead", "2g") \
    .config("spark.executor.memoryOverhead", "4g") \
    .config("spark.driver.cores", "4") \
    .config("spark.executor.cores", "4") \
    .config("spark.network.timeout", "1200s") \
    .config("spark.executor.heartbeatInterval", "120s") \
    .config("spark.sql.shuffle.partitions", "100") \
    .config("spark.default.parallelism", "8") \
    .config("spark.driver.maxResultSize", "0") \
    .getOrCreate()


    def _end_spark(self):
        """ Close the Spark session when all Spark operations are done """
        self.spark.stop()

    def _size_info(self, df):
        """ Get size information """
        sizes = df.rdd.map(lambda x: len(str(x))).collect()
        print(f"\n\nAverage size of rows: {sum(sizes) / len(sizes)}")
        print(f"Number of rows: {len(sizes)}")

    def _convert_lst_rdd(self, lst):
        """ Convert the list to an RDD """
        rdd = self.spark.sparkContext.parallelize(lst)
        return rdd

    def _load_from_parquet(self, parquet_dir):
        """ Read from parquet """
        if os.path.isdir(parquet_dir):
            parquet_files = glob.glob(os.path.join(parquet_dir, "*.parquet"))

            # Check if there are any .parquet files in the directory
            if parquet_files:
                print(f"\n\n{parquet_dir} is a folder and contains {len(parquet_files)} .parquet files.")
                return self.spark.read.parquet(parquet_dir)
            else:
                print(f"\n\n{parquet_dir} is a folder, but it contains no .parquet files.")
                raise ValueError("No .parquet files found")
        else:
            print(f"\n\n{parquet_dir} is not a folder.")
            raise ValueError("No .parquet files found")

    def _save_in_parquet(self, df_spark, path, _save=False):
        """
        Save a Spark DataFrame as a Parquet file with considerations for size and performance.

        :param df_spark: Spark DataFrame to be saved
        :param path: Path where the Parquet file will be saved
        :param _save: Boolean flag indicating whether to save the DataFrame
        """
        if not _save:
            return

        if not os.path.exists(path):
            os.makedirs(path)
        df_spark.write.mode("overwrite").parquet(path)

    def _damerau_levenshtein_distance(self, aligned_sequences):
        """
        This function implements Damerau-Levenshtein distance.

        Damerau-Levenshtein distance represents how many edits are needed to transform one sequence into another,
        where an edit can be an insertion, deletion, substitution, or transposition of adjacent characters (nucleotides).

        :param aligned_sequences:  Aligned sequences.
        :return: Pairwise Levenshtein distance.
        """

        seq_df = self.spark.createDataFrame([(seq,) for seq in aligned_sequences], ["sequence"])
        # seq_df.select("sequence").show(truncate=False)

        # Create a cross join to get all pairwise combinations
        cross_df = seq_df.crossJoin(seq_df.withColumnRenamed("sequence", "sequence2"))

        # For now, using a standard UDF as per your code
        distances_df = cross_df.withColumn("distance",
                                           F.udf(lambda s1, s2: textdistance.damerau_levenshtein(s1, s2))("sequence",
                                                                                                          "sequence2"))

        unique_distances_df = distances_df.filter(cross_df["sequence"] < cross_df["sequence2"])
        # unique_distances_df.select("distance").show(truncate=False)

        return unique_distances_df

    def _subgroup_dist_stats(self, distances_data):
        """
        For a subgroup, compute unique DNA barcode distance descriptive statistics such as min, max median, variance and standard deviation.

        NOTE that in biological contexts (e.g., DNA sequences), the variance of distances (e.g., pairwise genetic distances)
        can indicate genetic diversity within a group:
        - Low Variance in Genetic Distance: This implies that the sequences are quite similar to each other,
                                            indicating low genetic diversity.
        - High Variance in Genetic Distance: This suggests more genetic diversity, with some sequences being
                                             significantly different from others.

        :param distances_data:
        :return: Statistics dictionary.
        """

        # Check if repartitioning is necessary
        if distances_data.rdd.getNumPartitions() != self.num_partitions:
            distances_data = distances_data.repartition(self.num_partitions)

        # Cache the DataFrame if needed
        # distances_data.cache()

        # Calculate statistics
        results = distances_data.agg(
            F.mean("distance").alias("mean"),
            F.variance("distance").alias("variance"),
            F.stddev("distance").alias("stddev"),
            F.min("distance").alias("min"),
            F.max("distance").alias("max")
        ).collect()[0]

        # Uncache the DataFrame after use
        # distances_data.unpersist()

        return {
            "mean": round(results["mean"], 4),
            "variance": round(results["variance"], 4),
            "stddev": round(results["stddev"], 4),
            "min": results["min"],
            "max": results["max"],
        }

    def _subgroup_dists(self, subgroup, max_barcodes=1000, _subgroups_stats=False):
        """
        Compute pairwise distances of a subgroup.
        :param subgroup: The corresponding subgroup of a taxonomic rank.
        :param max_barcodes: Maximum number of barcodes randomly sampled in the subgroup.
        :param _subgroups_stats: If compute statistics on the fly. NOTE it makes processing slow; better skip it.
        :return: Name and pairwise distances of the subgroup.
        """

        name, barcodes = subgroup

        # Random sampling
        if (max_barcodes > 0) and (len(barcodes) > max_barcodes):
            sequences = random.sample(barcodes, max_barcodes)
        else:  # species
            sequences = barcodes

        # Perform alignment
        aligned_sequences = self.mafft_align.perform_mafft_alignment(sequences, name)

        # Compute pairwise distances
        distances = self._damerau_levenshtein_distance(aligned_sequences)

        # Calculate the statistics
        if _subgroups_stats:
            stats = self._subgroup_dist_stats(distances)
            stats.update({"subgroup_name": name})
            dist_count = distances.count()
            stats.update({"Pairwise Distance Count": dist_count})

        return (name, distances)

    def _rank_dists(self, rank_hierarchy, rank, min_barcodes=4, max_barcodes=1000, chunk_size=1000, path=None):
        """
          This function process distance computation per taxonomic rank.
          :param rank_hierarchy: Data hierarchy at a specified taxonomic level.
          :param rank: Taxonomic group level (e.g., family, genus, species).
          :param min_barcodes: Minimum number of barcodes per subgroups to compute pairwise distances.
          :param max_barcodes: Maximum number of barcodes per subgroups randomly sampled to compute pairwise distances.
          :param chunk_size: Chunk size of the subgroups of the rank.
          """

        tuple_list = [
            (subgroup, list(subgroup_entry['barcodes'].keys()))
            for subgroup, subgroup_entry in rank_hierarchy.items()
            if len(subgroup_entry['barcodes']) > min_barcodes
        ]

        tuple_chunks = [tuple_list[i:i + chunk_size] for i in range(0, len(tuple_list), chunk_size) ]

        chk_ended = 0
        max_barcodes = 0 if rank == 'species' else max_barcodes  # Process all barcodes of the species
        for chk in range(len(tuple_chunks)):

            if chk < chk_ended:
                continue

            chunk = tuple_chunks[chk]

            final_distances = None
            for cnt, item in tqdm(enumerate(chunk), total=len(chunk), desc=f"Processing subgroups of {rank}"):

                name, spark_df = self._subgroup_dists(item, max_barcodes=max_barcodes)

                # ----- If distances as spark_df
                spark_df = spark_df.withColumn("subgroup_name", F.lit(name))
                distances = spark_df.select("distance", "subgroup_name")
                distances = distances.withColumn("distance", col("distance").cast("float"))
                final_distances = distances if final_distances is None else final_distances.union(distances)

            # ---- Save the spark dataframe of pairwise distance for the subgroups existed in the chunk chk
            path = path if path is not None else self.save_path
            distances_path = os.path.join(path, rank, f"chunk_{chk}")
            if not os.path.exists(distances_path):
                os.makedirs(distances_path)
            self._save_in_parquet(final_distances, distances_path, _save=True)
            print(f'\n{rank} pairwise distances saved to {distances_path}.\n')

    def _rank_dist_stats(self, rank, chunks, distances_root=None, save_distances_pandas=False):
        """
        Compute pairwise distance statistics across taxonomic levels.
        :param rank: Taxonomic group level (e.g., family, genus, species).
        :param chunks: Available chunks of the subgroups of the rank.
        :param distances_root: Root path of the distances directory.
        :param save_distances_pandas: Covert and save to Pandas dataframe as a parquet file.
        """

        rank_stats = None
        df_spark = None
        for chk_num, chk in tqdm(enumerate(chunks), total=len(chunks), desc=f"Processing statistics of {rank}"):

            distances_dir = os.path.join(distances_root, f"chunk_{chk}")

            df = self._load_from_parquet(distances_dir)
            df_spark = df if df_spark is None else df_spark.union(df)

            # Group by 'subgroup_name' and calculate statistics for each group of the chunk
            stats_df = df.groupBy("subgroup_name").agg(
                F.mean("distance").alias("mean"),
                F.variance("distance").alias("variance"),
                F.stddev("distance").alias("stddev"),
                F.min("distance").alias("min"),
                F.max("distance").alias("max")
            )
            rank_stats = stats_df if rank_stats is None else rank_stats.union(stats_df)

        # Now aggregate the statistics (mean aggregation function) across the entire rank_stats DataFrame
        aggregated_rank_stat = rank_stats.agg(
            F.mean("mean").alias("mean_of_means"),
            F.mean("variance").alias("mean_of_variance"),
            F.mean("stddev").alias("mean_of_stddev"),
            F.mean("min").alias("mean_of_min"),
            F.mean("max").alias("mean_of_max")
        )

        # Show the result
        # aggregated_rank_stat.show(truncate=False)

        # Collect the results into a row
        aggregated_row = aggregated_rank_stat.collect()[0]

        # Convert the row to a dictionary
        rank_stats_dict = {
            "PWD-Mean": aggregated_row["mean_of_means"],
            "PWD-Variance": aggregated_row["mean_of_variance"],
            "PWD-Std.Dev": aggregated_row["mean_of_stddev"],
            "PWD-Min": aggregated_row["mean_of_min"],
            "PWD-Max": aggregated_row["mean_of_max"]
        }

        # Additionally convert and save distances in Pandas dataframe only if necessary: timely expensive
        path_pd = os.path.join(self.save_path, f"barcodes_pwd_{rank}.csv")
        save_in_pandas(df_spark, path_pd, _save=save_distances_pandas)

        return rank_stats_dict


