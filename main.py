
import os
import argparse
from barcode_statistic import BarcodeMetric


def main(config):

    barmetr = BarcodeMetric(metadata_file=config.metadata_file, method=config.method)

    # Create data hierarchy from metadata
    data_hierarchy = barmetr.build_hierarchy()

    # Compute Shannon Diversity Index (SDI) inplace
    data_hierarchy = barmetr.compute_sdi(data_hierarchy)

    # Save in pickle to prevent re-computation
    barmetr.create_pickle(data_hierarchy, os.path.join(barmetr.save_path, "data_hierarchy.pkl"))

    # Compute and save in parquets identical DNA pairwise distances
    barmetr.compute_pwd(data_hierarchy)

    # Compute the full statistics of the identical DNA barcodes
    barcode_stats = barmetr.compute_full_statistics(data_hierarchy)

    # Print full statistics
    barmetr.print_table(barcode_stats,
                        title="Full Statistics of Identical DNA Barcodes",
                        display_table=True)

    barmetr.save_in_pandas(barcode_stats, os.path.join(barmetr.save_path, "barcode_stats.csv"))



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compute DNA Barcode Statistics")

    parser.add_argument("--metadata_file", type=str, help="Path to the metadata CSV file",
                        default="/home/zahra/Desktop/PhyloTransNet/tmp/bioscan5m/metadata/csv/BIOSCAN_5M_Insect_Dataset_metadata.csv")
    parser.add_argument("--method", type=str, choices=["pandas", "spark"], default="spark",
                        help="Method to use for processing (pandas or spark)")

    # Parse the arguments
    args = parser.parse_args()

    main(args)

