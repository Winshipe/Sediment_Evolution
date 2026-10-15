# Data

Here's all the files you need to repeat my analyses. Large text files are gzipped; R's
`read_delim()` reads `.gz` files directly, so there is no need to unpack them.

## Naming conventions

All files use the same naming convention more or less so here's a breakdown of the short hand:

| Short hand | Meaning |
|---|---|
| `IS` | InStrain (the SNP & analysis pipeline, https://instrain.readthedocs.io/en/latest/) |
| `M1` / `M5` | the sample sites M1 or M5 |
| `M5_3`, `M5_25` | sample site plus MAG origin (MAGs were generated at 3 cm below sea floor (cmbsf) and 25 cmbsf) |
| `M5_3_10` | sample site, MAG origin, and read origin depth (in cmbsf) |
| `raw` | not normalised |
| `proc` | filtered (according to IS standards) and normalised (w/ https://github.com/Winshipe/Normalise-Paired-Regions) |
| `stb` | scaffold to bin, associates contig names w/ the MAGs |

"depth" always means sediment depth (cmbsf); depth of coverage is called "coverage".

## Folders

### `intermediate/` – tables passed between analysis steps

| File | Made by | Used by |
|---|---|---|
| `raw_genes_filtered_w_theta_hat_M5_3.tsv.gz`, `..._M1_3.tsv.gz` | InStrain gene_info for every read depth, plus the number of synonymous (`Ls`) and non-synonymous (`Ln`) sites per gene | `01` |
| `MAG_popgen_plus_anno_M5.tsv.gz`, `..._M1.tsv.gz` | `01` | `03` |
| `MAG_table_2.tsv` | `01` | `03` |

Running `01_make_gene_dataframes.Rmd` overwrites the `MAG_*` files here.

### `all_depths_M5_3/`, `all_depths_M1_3/` – InStrain output

`IS_<raw|proc>_<MAG set>_<read depth>_<SNVs|gene_info|genome_info>.tsv.gz`, one per read depth,
from `pipeline/Snakefile`. `all_depths_M5_3/` also contains:

* `combined_reference_*.faa.gz` / `.stb.gz` – predicted proteins and scaffold-to-bin file of the combined reference
* `clustered_aaai50_cluster.tsv.gz` / `_rep_seq.fasta.gz` – M5 genes clustered w/ `mmseqs easy-cluster` at 50% amino acid identity (if you've got something that really works to cluster these things in a meaningful way please message me on linkedin!)
* `file2family.txt.gz` – limited taxonomic information from GTDB
* `combined_annotations.txt.gz` – eggNOG-mapper output for the M5_3 MAGs
* `seq_for_analysis.fa.gz` – protein sequences

### `annotations/` – gene annotations

| File | Content |
|---|---|
| `uniref_annotations_*.tsv.gz` | UniProt ID-mapping output for the best UniRef90 hit of each gene |
| `phrog_hits_*.tsv.gz` | MMseqs2 hits against PHROGs (phage protein DB), BLAST tab format, no header |
| `phrog_table.tsv.gz` | PHROG annotation table |
| `mobileOG_annotations_M1_M5.tsv.gz` | hits against mobileOG-db (mobile genetic element proteins) |
| `busco_gene_ids_*.txt` | genes identified as BUSCO (single-copy core) genes |
| `eggnog_combined_annotations.txt` | eggNOG-mapper (https://github.com/eggnogdb/eggnog-mapper) output |

**Not included:** `uniref_hits_M5_3.tsv` and `uniref_hits_M1_3.tsv`, the raw MMseqs2 UniRef90
hit tables needed to rerun `01_make_gene_dataframes.Rmd` (made by the `mmseqs_uniref` rule
in `pipeline/Snakefile`).

### `mags/` – MAG metadata

| File | Content |
|---|---|
| `MAG_table_further_analysis.tsv` | taxonomy, BUSCO and CheckM2 statistics of the focal MAGs |
| `checkm2_quality_report_*.tsv` | CheckM2 report for the focal MAGs |
| `checkm2_results_all_*.tsv` | CheckM2 report for all MAGs |
| `*.gtdbtk.bac120.summary.tsv` | GTDB-Tk classification of all MAGs |
| `lifestyles_*.tsv` | lifestyle (surface specialist / generalist / subsurface specialist) and short phylogeny label of each focal MAG |
| `combined_sites_bacteria_tree_vft.nwk` | phylogenomic tree of all bacterial MAGs from both sites |

### `coverages/`

`samtools coverage` output per scaffold for every MAG set / read depth, plus the scaffold-to-bin
files to assign scaffolds to MAGs. Used for the tree in Figure 1.

### `amplicon_16S/`

16S rRNA amplicon data used in Figure 1 (B, C): `1ASV_raw_abund_tax.txt` (ASV table with
taxonomy) and `Lilian_metadata.txt` (sample metadata).

### `images/`

`Sediment_core_cartoon_v2.svg` – the sediment core drawing used as Figure 1A.

### `other/`

Older InStrain runs and an older tree, not used by the current notebooks.
