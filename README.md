# barcodemetr
This repository contains code for analyzing DNA nucleotide barcode sequences. 
The DNA nucleotide sequences are organized by taxonomic ranks—phylum, class, order, family, subfamily, genus, and species. 
At each rank, the data is structured hierarchically: every taxonomic rank consists of subgroups, and each subgroup 
includes identical DNA barcode sequences. Figure 1 illustrates this hierarchical organization.

### DNA Barcode Sequence

The DNA barcode sequence highlights the arrangement of the four nucleotides—Adenine (A), Thymine (T), Cytosine (C), 
and Guanine (G)—within a specific gene region, such as the mitochondrial cytochrome c oxidase subunit I (COI) gene. 
Each nucleotide is depicted using a distinct color block to enhance visual interpretation:

```
TTTATATTTTATTTTTGGAGCATGATCAGGAATAGTTGGAACTTCAATAAGTTTATTAATTCGAACAGAATTAAGCCAACCAGGATCAACATTTAT ....
```

- Adenine (A): Red
- Thymine (T): Blue
- Cytosine (C): Green
- Guanine (G): Yellow


### Damerau-Levenshtein Distance for DNA Barcodes
The Damerau-Levenshtein distance extends the standard Levenshtein metric by including adjacent character transpositions 
in addition to insertions, deletions, and substitutions. It quantifies the similarity between DNA sequences based on 
the minimum number of such edits required to transform one sequence into another.

The code computes the average pairwise distances between unique DNA barcodes at various taxonomic ranks. For each subgroup:

- Groups with <4 unique sequences are excluded.
- If >1,000 sequences exist, 1,000 are randomly sampled.
- Sequences are first aligned using MAFFT.
- Pairwise distances are computed using the Damerau-Levenshtein metric.
- Mean and standard deviation are aggregated across subgroups.

### 
The repository supports two implementations of DNA barcodes pairwise distance calculations one with Pandas and the other with Apache Spark (PySpark).
For optimal performnace conduct experiments in two separate phases:

> **ℹ️ Note** Use the Pandas implementation for small datasets. For large datasets such as BIOSCAN-5M, the Apache Spark implementation is recommended.

#### I. Calculate&Save SDI and Pairwise Distances
To initially calculate DNA barcodes pairwise distances and save them in the Parquet files or Pandas dataframes, execute:

```bash
python main.py --method spark --compute_pwd --load_metadata --metadata_file <file-path> --ranked_data_file <file-path> --save_path <directory-path>
``` 

#### II.  Enable Statistical Processing of SDI and Pairwise Distances
To extract the statistics based on saved sdi and pairwise distances execute:

```bash
python main.py --method spark --compute_full_statistics --display_table --save_statistics --ranked_data_file <file-path> --save_path <directory-path>
``` 

#### III.  Visualization
To visualize the statics execute the following:

```bash
python visualization.py --create_plots --rank <Taxonomic-rank> --metadata_file <file-path> --ranked_data_file <file-path> --distances_path <directory-path> --save_path <directory-path>
``` 
<div align="center">
  <figure>
    <img src="figures/sdi_distributions.png" 
         alt="class." />
    <figcaption><b>Figure 1: </b>Distribution of Shannon Diversity Index (SDI) across subgroups of taxonomic ranks. 
          The x-axis shows the Taxonomic ranks.
    categories sorted alphabetically.</figcaption>
  </figure>
</div>

<div align="center">
  <figure>
    <img src="figures/class_distance_distribution.png" 
         alt="class." />
    <figcaption><b>Figure 1: </b>Distribution of pairwise distances of subgroups of class. The x-axis shows the subgroup
    categories sorted alphabetically.</figcaption>
  </figure>
</div>

<div align="center">
  <figure>
    <img src="figures/order_distance_distribution.png" 
         alt="order." />
    <figcaption><b>Figure 2: </b>Distribution of pairwise distances of subgroups of order. The x-axis shows the subgroup
    categories sorted alphabetically.</figcaption>
  </figure>
</div>

<div align="center">
  <figure>
    <img src="figures/species_distance_distribution.png" 
         alt="species." />
  </figure>
</div>
<div align="left">
  <p><b>Figure 3:</b> Distribution of pairwise distances of subgroups of species. Among the species, there are
    8,372 distinct subgroups with sufficient identical barcodes for calculating pairwise distances, which
    makes visualization challenging. To address this, the groups are sorted in descending order based
    on their mean distances and partitioned into 100 bins. These bins are used to plot the distribution
    of pairwise distances within the species rank. The mean distance of each bin is displayed along the
    x-axis.</div>


