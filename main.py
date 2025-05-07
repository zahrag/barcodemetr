
import os
import argparse
from barcode_statistic import BarcodeMetric
from utils import *

def main(configs):

    barmetr = BarcodeMetric(metadata_file= configs.metadata_file,
                            method=configs.method, load_metadata=configs.load_metadata)

    # Create ranked data hierarchy from metadata
    ranked_data = barmetr.build_hierarchy(path=configs.ranked_data_path)

    # Compute Shannon Diversity Index (SDI) inplace
    ranked_data = barmetr.compute_sdi(ranked_data, enabled=configs.compute_sdi)

    # Save in pickle to prevent re-computation
    create_pickle(data=None, pickle_file="")

    # Compute and save in parquets identical DNA pairwise distances
    barmetr.compute_pwd(ranked_data=None)

    # Compute the full statistics of the identical DNA barcodes
    barcode_stats_full = barmetr.compute_full_statistics(ranked_data, save_distances_pandas=configs.save_distances_pandas)

    # Print full statistics
    print_table(barcode_stats_full, title="Full Statistics of Identical DNA Barcodes", display_table=configs.display_table
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
                        default="/home/zahra/Desktop/PhyloTransNet/tmp/bioscan5m/metadata/csv/BIOSCAN_5M_Insect_Dataset_metadata.csv")

    parser.add_argument("--method",
                        type=str,
                        choices=["pandas", "spark"],
                        default="spark",
                        help="Method to use for processing (pandas or spark)")

    parser.add_argument("--ranked_data_path",
                        type=str,
                        default="/home/zahra/Desktop/PhyloTransNet/tmp/bioscan5m/metadata/csv/barcode_analysis/spark/data_hierarchy.pkl",
                        help="Path to the ranked data pickle file",)

    parser.add_argument('--load_metadata',
                        default=False,
                        action='store_true',
                        help='IF loading metadata (save time by setting False when ranked data is available)?')


    parser.add_argument('--compute_sdi',
                        default=False,
                        action='store_true',
                        help='IF computing Shannon Diversity Index (SDI)?')

    parser.add_argument('--compute_pwd',
                        default=False,
                        action='store_true',
                        help='IF computing pairwise distances of identical DNA barcode sequences?')

    parser.add_argument('--display_table',
                        default=False,
                        action='store_true',
                        help='IF displaying table of full DNA barcode statistics?')

    parser.add_argument('--save_statistics',
                        default=False,
                        action='store_true',
                        help='IF saving the full DNA barcode statistics?')

    parser.add_argument('--compute_full_statistics',
                        default=False,
                        action='store_true',
                        help='IF saving the pairwise distance of each rank in panda? (NOTE that it is time consuming.)')

    # Parse the arguments
    args = parser.parse_args()

    main(args)

