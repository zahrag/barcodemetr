
import matplotlib.pyplot as plt
import os


def plot_box(df, fig_file,
             color_dict=None,
             custom_order=None,
             fig_format='png',
             plt_type = 'box',
             x="subgroup_name",
             y="distance",
             x_label='Group Name',
             y_label='Distance',
             tick_font_size=24,
             _plt=False):

    if not _plt:
        return

    import plotly.express as px
    import time

    fig = px.box(df, x=x, y=y, color=x, color_discrete_map=color_dict)

    fig.update_layout(
        xaxis=dict(
            categoryorder='array',
            categoryarray=custom_order,
            tickangle=90,  # Rotate x-axis labels by 90 degrees
            title=dict(
                text=x_label,  # Set the x-axis label
                font=dict(size=tick_font_size)  # Set font size for x-axis label
            ),
            tickfont=dict(size=tick_font_size)
        ),
        yaxis=dict(
            title=dict(
                text=y_label,  # Set the y-axis label
                font=dict(size=tick_font_size)  # Set font size for y-axis label
            ),
            tickfont = dict(size=tick_font_size)
        ),
        showlegend=False
    )
    plt.tight_layout()

    time.sleep(1)

    fig_path = os.path.join(os.path.dirname(fig_file), plt_type)
    fig_name = os.path.splitext(os.path.basename(fig_file))[0]

    fig.show()
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)

    # Save the Plotly figure as a PDF
    if fig_format == 'pdf':
        fig.write_image(os.path.join(fig_path, f'{fig_name}_{plt_type}.pdf'),
                        format="pdf",
                        width=1920,  # Width in pixels
                        height=1080,  # Height in pixels
                        scale=1)

    elif fig_format == 'png':
        fig.write_image(os.path.join(fig_path, f'{fig_name}_{plt_type}.png'),
                        format="png",
                        width=1920,
                        height=1080,
                        scale=3)
    elif fig_format == 'svg':
        fig.write_image(os.path.join(fig_path, f'{fig_name}_{plt_type}.svg'),
                        format="svg",
                        width=1920,
                        height=1080,
                        scale=1)

    fig.show()


def plot_box_bin(df, fig_file,
                 fig_format='pdf',
                 plt_type = 'box',
                 x="bins",
                 y="mean_distance",
                 x_label='Bins',
                 y_label='Mean Distance',
                 tick_font_size=16,
                 _plt=False,
                 ):

    if not _plt:
        return

    import plotly.express as px

    # Create the box plot
    fig = px.box(df, x=x, y=y)

    # Calculate mean distances for each bin
    bin_means = df.groupby(x)[y].mean().reset_index()

    # Extract mean distances
    mean_distances = bin_means[y].values

    # Extract bin labels
    bin_labels = bin_means[x].values

    # Update x-axis ticks and labels
    fig.update_layout(
        xaxis=dict(
            tickvals=bin_labels,  # Set tick values to the bins
            ticktext=[f"{label}: {mean:.2f}" for label, mean in zip(bin_labels, mean_distances)],
            tickfont=dict(size=tick_font_size)
        ),
        yaxis=dict(
            tickfont=dict(size=tick_font_size)
        ),

    )

    # Update layout
    fig.update_layout(
        xaxis=dict(
            categoryorder='total descending',
            tickangle=90,
            title=dict(
                text=x_label,
                font=dict(size=24)
            )
        ),
        yaxis=dict(
            title=dict(
                text=y_label,
                font=dict(size=24)
            )
        ),
        showlegend=False
    )
    plt.tight_layout()

    # Define output paths
    fig_path = os.path.join(os.path.dirname(fig_file), plt_type)
    fig_name = os.path.splitext(os.path.basename(fig_file))[0]

    if fig_format == 'pdf':
        fig.write_image(os.path.join(fig_path, f'{fig_name}_{plt_type}.pdf'),
                        format="pdf",
                        width=1920,  # Width in pixels
                        height=1080,  # Height in pixels
                        scale=1)

    elif fig_format == 'png':
        fig.write_image(os.path.join(fig_path, f'{fig_name}_{plt_type}.png'),
                        format="png",
                        width=1920,
                        height=1080,
                        scale=3)
    elif fig_format == 'svg':
        fig.write_image(os.path.join(fig_path, f'{fig_name}_{plt_type}.svg'),
                        format="svg",
                        width=1920,
                        height=1080,
                        scale=1)
    fig.show()

