# barcodemetr
This repository contains code for analyzing DNA nucleotide barcode sequences. 
The DNA nucleotide sequences are organized by taxonomic ranks—phylum, class, order, family, subfamily, genus, and species. 
At each rank, the data is structured hierarchically: every group-level consists of subgroups, and each subgroup 
includes identical DNA barcode sequences. Figure 1 illustrates this hierarchical organization.

### DNA Barcode Sequence

The illustrated DNA barcode sequence highlights the arrangement of the four nucleotides—Adenine (A), Thymine (T), Cytosine (C), 
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


