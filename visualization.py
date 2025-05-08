
import argparse
from utils import *
from plot_functions import *
from barcode_pwd_spark import BarcodePWD
from tqdm import tqdm
import matplotlib.colors as mcolors


def _plt_sdi(args, ranked_data, _plt=False):

    if not _plt:
        return

    # Example: flatten data_hierarchy into a DataFrame
    rows = []
    for rank, subgroups in ranked_data.items():
        for subgroup, info in subgroups.items():
            rows.append({'rank': rank, 'subgroup': subgroup, 'sdi': info['sdi']})

    df_sdi = pd.DataFrame(rows)

    weights = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
    colors = plt.cm.Reds(weights)
    color_hex = [mcolors.to_hex(color) for color in colors]
    ranks = ["phylum", "class", "order", "family", "subfamily", "genus", "species"]
    color_dict = dict(zip(ranks, color_hex))

    custom_order = list(color_dict.keys())

    fig_file = os.path.join(args.save_path, f'sdi_distributions.pdf')
    plot_box(df_sdi, fig_file,
             color_dict=color_dict,
             custom_order=custom_order,
             x="rank",
             y="sdi",
             fig_format='png',
             x_label='Taxonomic Ranks',
             y_label='Shannon Diversity Index',
             tick_font_size=24,
             _plt=_plt)

def load_distances(args, _to_panda=False):

    pwd = BarcodePWD()
    distances_dir = os.path.join(args.distances_path, args.rank)
    chks = extract_chunks(args.rank, distances_dir, method="spark")
    print(f"\nSaved chunks: {chks}.\n")
    distances = None
    for chk_num, chk in tqdm(enumerate(chks), total=len(chks), desc=f"Processing statistics of {args.rank}"):
        distances_dir_chk = os.path.join(distances_dir, f"chunk_{chk}")
        df = pwd._load_from_parquet(distances_dir_chk)
        distances = df if distances is None else distances.union(df)

    # Convert Spark DataFrame to Pandas: NOTE that this takes time
    distances_pandas = None
    if _to_panda:
        distances_pandas = distances.toPandas()

    return distances, distances_pandas

def check_match_mismatch(args, distances, ranked_data, _check=False):

    if not _check:
        return

    for rank, subgroups in ranked_data.items():
        if rank != args.rank:
            continue
        for subgroup_name, subgroup_data in subgroups.items():
            n_barcodes = len(subgroup_data['barcodes'])
            if n_barcodes < 5:
                continue
            # Expected number of comparisons
            n = 1000 if len(subgroup_data['barcodes']) > 1000 else len(subgroup_data['barcodes'])
            expected = n * (n - 1) // 2

            # Actual number of comparisons in your DataFrame
            if isinstance(distances, pd.DataFrame):
                actual = distances[distances['subgroup_name'] == subgroup_name].shape[0]
            else:
                actual = distances.filter(distances['subgroup_name'] == subgroup_name).count()

            print(f"[{rank}::{subgroup_name}] expected={expected}, actual={actual}")
            if expected == actual:
                print("  ✅ Match")
            else:
                print("  ❌ Mismatch")


def agg_bin_subgroups(df, num_bins=100, _plt=False):

    # Get a df where each entry is a unique group of df (under 'group_name'), and its mean distance (under 'distance')
    mean_distances = df.groupby('subgroup_name')['distance'].mean().reset_index()

    # Rename the columns for clarity: 'distance' changed to 'mean_distance'
    mean_distances.columns = ['subgroup_name', 'mean_distance']

    # Rank the groups based on mean distance descending (Larger numbers come first, followed by smaller numbers)
    mean_distances['rank'] = mean_distances['mean_distance'].rank(method='first', ascending=False)

    # Sort by rank for better readability
    mean_distances = mean_distances.sort_values(by='rank')

    mean_distances['bins'] = pd.qcut(mean_distances['rank'],
                                     q=num_bins,
                                     labels=[f'Bin {i + 1}' for i in range(num_bins)]
                                     )

    # Calculate mean distances for each bin
    bin_means = mean_distances.groupby('bins')['mean_distance'].mean()

    if _plt:
        plt.figure(figsize=(10, 6))
        mean_distances.groupby('bins')['mean_distance'].mean().plot(kind='bar')  # Fix here
        plt.title('Mean Distances by Bins')
        plt.xlabel('Bins')
        plt.ylabel('Mean Distance')
        # Set x-ticks to show the mean distances instead of Bin 1, Bin 2, ...
        plt.xticks(ticks=range(len(bin_means)), labels=[f'{mean:.2f}' for mean in bin_means], rotation=90, fontsize=7)
        plt.tight_layout()
        plt.show()

    return mean_distances


def main(args):

    ranked_data = open_pickle(pickle_file=args.ranked_data_file)

    # Plot SDI distributions Across 7 Taxonomic Ranks
    _plt_sdi(args, ranked_data, _plt=False)

    # Load pairwise distances of DNA barcode from parquets
    distances_spark, distances = load_distances(args, _to_panda=True)

    # Sanity checks of the number of distances
    check_match_mismatch(args, distances_spark, ranked_data, _check=False)

    if args.rank in ['family', 'subfamily', 'genus', 'species']:

        assert isinstance(distances, pd.DataFrame), "Only Pandas DataFrames are supported for plotting."

        unique_groups = distances['subgroup_name'].unique()
        print(f'Number of unique subgroups of {args.rank}: {len(unique_groups)}')

        num_bins = 10 #100
        df = agg_bin_subgroups(distances, num_bins=num_bins)

        fig_file = os.path.join(args.save_path, f'{args.rank}_distributions.pdf')
        plot_box_bin(df, fig_file,
                     fig_format='png',
                     x="bins",
                     y="mean_distance",
                     x_label=f'{args.rank.capitalize()} name',
                     y_label='Damerau-Levenshtein',
                     tick_font_size=18,
                     _plt=args.create_plots)

    elif args.rank in ["class", "order"]:

        assert isinstance(distances, pd.DataFrame), "Only Pandas DataFrames are supported for plotting."

        unique_groups = distances['subgroup_name'].unique()
        print(f'Number of unique subgroups of {args.rank}: {len(unique_groups)}')

        rank_dict = ranked_data[args.rank]
        # Sort the subgroup names within that rank by descending SDI value.
        sorted_group_names = sorted(rank_dict.keys(), key=lambda g: rank_dict[g]['sdi'], reverse=True)
        # Re-sorts the same group names alphabetically.
        sorted_groups_alphabetically = sorted(sorted_group_names)  # or names.sort(reverse=True)

        # Load metadata file
        metadata = load_from_pandas(args.metadata_file, load_file=True)

        # Define a list of 10 colors
        color_palette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                         "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",]

        # Get unique classes from the 'class' column
        unique_classes = metadata['class'].dropna().unique()

        # Create a dictionary mapping each class to a unique color
        class_color_map = {class_name: color_palette[i] for i, class_name in enumerate(unique_classes)}

        # Filter the DataFrame to exclude rows with NaN in 'class' or 'order'
        filtered_df = metadata.dropna(subset=['class', 'order'])

        # Map colors to each order based on its class
        order_color_map = {order: class_color_map[class_name]
                           for order, class_name in zip(filtered_df['order'], filtered_df['class'])}

        # Map order to their corresponding classes
        order_class_map = {order: class_name for order, class_name in zip(filtered_df['order'], filtered_df['class'])}
        print(order_class_map)

        if args.rank == 'class':
            color_dict = class_color_map
        else:
            color_dict = order_color_map

        fig_file = os.path.join(args.save_path, f'{args.rank}_distributions.pdf')
        plot_box(distances, fig_file,
                 color_dict=color_dict,
                 custom_order=sorted_groups_alphabetically,
                 x="subgroup_name",
                 y="distance",
                 fig_format='png',
                 x_label=f'{args.rank.capitalize()} name',
                 y_label='Damerau-Levenshtein',
                 tick_font_size=24,
                 _plt=args.create_plots)
    else:
        raise ValueError(f'Invalid rank {args.rank}')



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Compute DNA Barcode Statistics")

    parser.add_argument("--rank",
                        type=str,
                        help="taxonomic rank to plot distributions for",
                        default="family")

    parser.add_argument("--metadata_file",
                        type=str,
                        help="Path to the metadata CSV file",
                        default="./BIOSCAN_5M_Insect_Dataset_metadata.csv")

    parser.add_argument("--ranked_data_file",
                        type=str,
                        help="Path to the ranked data Pickle file",
                        default="./data_hierarchy.pkl")

    parser.add_argument("--distances_path",
                        type=str,
                        help="Path to the distances Parquet files",
                        default="./distances")

    parser.add_argument("--save_path",
                        type=str,
                        help="Path to the directory to save the figures",
                        default="./distances_figures")

    parser.add_argument('--create_plots',
                        default=False,
                        action='store_true',
                        help='IF plotting?')

    configs = parser.parse_args()
    main(configs)





