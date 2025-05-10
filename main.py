
import os
import argparse
from barcode_statistic import BarcodeMetric
from utils import *

def main(configs):

    barmetr = BarcodeMetric(method=configs.method,
                            max_barcode_length=configs.max_barcode_length,
                            metadata_file= configs.metadata_file,
                            load_metadata=configs.load_metadata,
                            save_path=configs.save_path,
                            )

    # Create ranked data hierarchy from metadata
    ranked_data = barmetr.build_hierarchy(path=configs.ranked_data_file)

    # Compute Shannon Diversity Index (SDI) inplace
    ranked_data = barmetr.compute_sdi(ranked_data, path=configs.ranked_data_file)

    # Save in pickle to prevent re-computation
    create_pickle(ranked_data, path=configs.ranked_data_file)

    # Compute and save in parquets the identical DNA's pairwise distances
    barmetr.compute_pwd(ranked_data, _enabled=configs.compute_pwd)

    # Compute the full statistics of the identical DNA barcodes
    barcode_stats_full = barmetr.compute_full_statistics(ranked_data,
                                                         save_distances_pandas=configs.save_distances_pandas,
                                                         enabled=configs.compute_full_statistics,
                                                         )

    # Print full statistics
    print_table(barcode_stats_full,
                title="Full Statistics of Identical DNA Barcodes",
                display_table=configs.display_table,
                )

    save_in_pandas(barcode_stats_full,
                   os.path.join(barmetr.save_path, "barcode_stats.csv"),
                   _save=configs.save_statistics
                   )



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compute DNA Barcode Statistics")

    parser.add_argument("--metadata_file",
                        type=str,
                        help="Path to the metadata CSV file",
                        default="")

    parser.add_argument("--method",
                        type=str,
                        choices=["pandas", "spark"],
                        default="spark",
                        help="Method to use for processing (pandas or spark)")

    parser.add_argument("--ranked_data_file",
                        type=str,
                        default="",
                        help="Path to the ranked data pickle file",)

    parser.add_argument("--save_path",
                        type=str,
                        default="",
                        help="Path to save the results", )
    parser.add_argument("--max_barcode_length",
                        type=int,
                        default=625,
                        help="Maximum barcode length applied to align barcodes.", )

    parser.add_argument('--load_metadata',
                        default=False,
                        action='store_true',
                        help='IF loading metadata (save time by setting False when ranked data is available)?')

    parser.add_argument('--compute_pwd',
                        default=False,
                        action='store_true',
                        help='IF computing pairwise distances of identical DNA barcode sequences?')

    parser.add_argument('--compute_full_statistics',
                        default=False,
                        action='store_true',
                        help='IF enabling calculations of full DNA barcode statistics?')

    parser.add_argument('--display_table',
                        default=False,
                        action='store_true',
                        help='IF displaying table of full DNA barcode statistics?')

    parser.add_argument('--save_statistics',
                        default=False,
                        action='store_true',
                        help='IF saving the full DNA barcode statistics?')

    parser.add_argument('--save_distances_pandas',
                        default=False,
                        action='store_true',
                        help='IF saving the pairwise distance of each rank in pandas? (NOTE that it is time consuming.)')

    args = parser.parse_args()
    main(args)

