

import textdistance
import os
import glob

from attr.validators import max_len
from gitdb.fun import chunk_size
from tqdm import tqdm
import logging
import random
from pyspark.sql import SparkSession
from pyspark.sql import functions as F
from pyspark.sql.functions import col
from functools import reduce
from pyspark.sql.functions import pandas_udf
from pyspark.sql.types import FloatType
import textdistance

import pandas as pd
from barcode_alignment import perform_mafft_alignment
from barcode_sampling import BarcodeSampler
from pathlib import Path
current_directory = Path(__file__).parent


# @pandas_udf(FloatType())
# def damerau_levenshtein_udf(s1: pd.Series, s2: pd.Series) -> pd.Series:
#     return s1.combine(s2, lambda x, y: textdistance.damerau_levenshtein(x, y))

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

    def __init__(self):

        # try
        self.spark = self.initialize_spark()
        # Determine the number of available cores
        self.num_cores = self.spark.sparkContext.defaultParallelism

        # Set the number of partitions to 2-4 times the number of cores
        self.num_partitions = max(2, min(4 * self.num_cores, 100))  # Limit to a reasonable max, e.g., 100

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


    def end_spark(self):
        """ Close the Spark session when all Spark operations are done """
        self.spark.stop()

    def convert_lst_rdd(self, lst):
        """ Convert the list to an RDD """
        rdd = self.spark.sparkContext.parallelize(lst)
        return rdd

    def read_df_spark(self, parquet_dir):
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

    def get_size_info(self, df):
        """ Get size information """
        sizes = df.rdd.map(lambda x: len(str(x))).collect()
        print(f"\n\nAverage size of rows: {sum(sizes) / len(sizes)}")
        print(f"Number of rows: {len(sizes)}")

    def save_as_parquet(self, df_spark, path, _save=False):
        """
        Save a Spark DataFrame as a Parquet file with considerations for size and performance.

        :param df_spark: Spark DataFrame to be saved
        :param path: Path where the Parquet file will be saved
        :param _save: Boolean flag indicating whether to save the DataFrame
        """
        if not _save:
            return

        # df_spark.cache()

        # Optionally repartition to optimize for saving
        # df_spark = df_spark.repartition(100)

        # df_spark.write.option("compression", "gzip").mode("overwrite").parquet(path)

        df_spark.write.mode("overwrite").parquet(path)

        # df_spark.unpersist()

    def _pwd(self, aligned_sequences):
        """
        This function implements Damerau-Levenshtein distance using PySpark.

        Damerau-Levenshtein distance represents how many edits are needed to transform one sequence into another,
        where an edit can be an insertion, deletion, substitution, or transposition of adjacent characters (nucleotides).

        :param aligned_sequences:  Aligned sequences.
        :return: Pairwise Levenshtein distance.
        """

        seq_df = self.spark.createDataFrame([(seq,) for seq in aligned_sequences], ["sequence"])
        # seq_df.select("sequence").show(truncate=False)

        # Create a cross join to get all pairwise combinations
        cross_df = seq_df.crossJoin(seq_df.withColumnRenamed("sequence", "sequence2"))

        # Using Pandas UDF for better performance, if applicable
        # distances_df = cross_df.withColumn("distance", pandas_udf("sequence", "sequence2"))

        # For now, using a standard UDF as per your code
        distances_df = cross_df.withColumn("distance",
                                           F.udf(lambda s1, s2: textdistance.damerau_levenshtein(s1, s2))("sequence",
                                                                                                          "sequence2"))

        # distances_df = cross_df.withColumn("distance", damerau_levenshtein_udf("sequence", "sequence2"))

        unique_distances_df = distances_df.filter(cross_df["sequence"] < cross_df["sequence2"])
        # unique_distances_df.select("distance").show(truncate=False)

        return unique_distances_df

    def get_statistics(self, distances_data):
        """
        In biological contexts (e.g., DNA sequences), the variance of distances (e.g., pairwise genetic distances)
        can indicate genetic diversity within a group.
        Low Variance in Genetic Distance: This implies that the sequences are quite similar to each other, indicating low genetic diversity.
        High Variance in Genetic Distance: This suggests more genetic diversity, with some sequences being significantly different from others.

        :param distances_data:
        :return:
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

    def _subgroup_distances(self, item, max_len=1000):

        name, sequences = item

        # Random sampling
        sequences = random.sample(sequences, max_len) if len(sequences) > max_len else sequences

        # Perform alignment
        aligned_sequences = perform_mafft_alignment(sequences, name)

        # Compute pairwise distances
        distances = self._pwd(aligned_sequences)

        # Calculate the statistics
        # stats = self.get_statistics(distances)
        # stats.update({"group_name": group_name})
        # dist_count = distances.count()
        # stats.update({"Pairwise Distance Count": dist_count})

        return (name, distances)

    def _rank_distances(self, data_dict, rank, process=False):

        if not process:
            return

        # Convert the dictionary to a list of tuples [(species, sequences), ...]
        print(f'Processing DNA barcodes pairwise distances across {rank} ...')

        min_seq_num = 3
        tuple_list = [(group_name, data['unique_barcodes'])
                      for group_name, data in data_dict.items()
                      if len(data['unique_barcodes']) > min_seq_num ]

        chunk_size = 1000
        tuple_list_chunks = [tuple_list[i:i + chunk_size] for i in range(0, len(tuple_list), chunk_size)]

        chk_ended = 0
        max_seq_length = 1000
        for chk in range(len(tuple_list_chunks)):

            if chk < chk_ended:
                continue

            chunk = tuple_list_chunks[chk]

            final_distances = None
            # stats_lst = []
            # distances_dict = {}
            for cnt, item in tqdm(enumerate(chunk), total=len(chunk), desc="Processing groups"):

                group_name, spark_df = self._subgroup_distances(item, max_len=max_seq_length)

                # ------ If statistics
                # stats_lst.append(stats)

                # ------ If distances as list
                # distances_dict[group_name] = spark_df.select("distance").rdd.flatMap(lambda x: x).collect()
                # distances_dict[group_name] = spark_df.select("distance").toPandas()["distance"].tolist()

                # ----- If distances as spark_df
                spark_df = spark_df.withColumn("group_name", F.lit(group_name))
                distances = spark_df.select("distance", "group_name")
                distances = distances.withColumn("distance", col("distance").cast("float"))
                final_distances = distances if final_distances is None else final_distances.union(distances)

            # ---- Save the DataFrame to Parquet format
            distances_dir = f"{current_directory}/distances/{rank}/chunk_{chk}"
            make_directory(distances_dir)
            # create_pickle(distances_dict, distances_dir)
            self.save_as_parquet(final_distances, distances_dir, _save=True)
            print('')


    def _rank_statistics(self, rank, max_chk=10, process=False):

        if not process:
            return None

        print(f'Processing DNA barcodes statistics across {rank} ...')

        rank_stats = None
        df_spark = None
        for chk in tqdm(range(max_chk), total=max_chk, desc="Processing statistics"):

            distances_dir = f"{current_directory}/distances/{rank}/chunk_{chk}"
            if not os.path.exists(distances_dir):
                continue

            df = self.read_df_spark(distances_dir)
            df_spark = df if df_spark is None else df_spark.union(df)

            # Group by 'group_name' and calculate statistics for each group of the chunk
            stats_df = df.groupBy("group_name").agg(
                F.mean("distance").alias("mean"),
                F.variance("distance").alias("variance"),
                F.stddev("distance").alias("stddev"),
                F.min("distance").alias("min"),
                F.max("distance").alias("max")
            )
            rank_stats = stats_df if rank_stats is None else rank_stats.union(stats_df)

        df_pandas = df_spark.toPandas()
        make_tsv(df_pandas, name=f'BIOSCAN_5M_distances_{rank}.tsv', path=os.path.join(current_directory, 'distances'))

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
            "Mean": aggregated_row["mean_of_means"],
            "Variance": aggregated_row["mean_of_variance"],
            "Std.Dev": aggregated_row["mean_of_stddev"],
            "Min": aggregated_row["mean_of_min"],
            "Max": aggregated_row["mean_of_max"]
        }

        return rank_stats_dict


