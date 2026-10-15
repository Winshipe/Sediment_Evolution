# Sediment Evolution

Population genetics of bacterial MAGs (metagenome-assembled genomes) along marine sediment
cores from two sites (M1, M5). Reads from sediment depths of 1–25 cm below sea floor (cmbsf) are
mapped back to MAGs assembled from 3 cmbsf, and per-gene πN/πS, θN/θS and pN/pS are compared
across depths to see how the efficacy of selection changes with depth.

## Layout

```
Sediment_Evolution/
├── analysis/                          R Markdown notebooks (run in this order)
│   ├── 01_make_gene_dataframes.Rmd    per-gene popgen stats + annotations -> data/intermediate/
│   ├── 02_figure_1_community.Rmd      Figure 1 (16S community, MAG tree, MAG coverage)
│   └── 03_figures_selection.Rmd       Figures 2-4, supplementary figure and tables
├── data/                              all inputs, see data/README.md
├── pipeline/                          Snakemake pipeline: reads -> MAG mapping -> InStrain (HPC)
│   ├── Snakefile
│   ├── config.yaml
│   └── envs/                          conda environments
├── results/                           created when the notebooks run (figures/, tables/)
├── install_R_packages.R
└── Sediment_Evolution.Rproj
```

All paths in the notebooks are relative to the project root (via the `here` package), so the
project works wherever it is unpacked.

## Quick start (figures)

Requires R ≥ 4.2. Knitting to HTML also needs pandoc, which comes with RStudio (from a plain
terminal, install pandoc separately).

1. Install the R packages (once):

   ```bash
   Rscript install_R_packages.R
   ```

2. Open `Sediment_Evolution.Rproj` in RStudio and knit / run the notebooks in `analysis/`,
   or from a terminal in the project folder:

   ```bash
   Rscript -e 'rmarkdown::render("analysis/03_figures_selection.Rmd")'
   ```

   ```bash
   Rscript -e 'rmarkdown::render("analysis/02_figure_1_community.Rmd")'
   ```

Figures are written to `results/figures/`.

`02` and `03` run directly from the shipped data. `01` regenerates the tables in
`data/intermediate/` that `03` reads; pre-computed copies are included, so you only need it
if you want to rebuild them. It needs two extra files (`data/annotations/uniref_hits_*.tsv`,
see `data/README.md`).

The bootstrap steps use random numbers, so values in the bootstrap tables differ slightly
between runs.

## Upstream pipeline (optional)

`pipeline/Snakefile` produces the InStrain output in `data/all_depths_*` from raw reads.
It needs the raw reads and MAG fasta files (not included, they are too large), plus
eggNOG-mapper and UniRef90 MMseqs2 databases and the
[Normalise-Paired-Regions](https://github.com/Winshipe/Normalise-Paired-Regions) tool. Set
their locations in `pipeline/config.yaml`, then:

```bash
cd pipeline && snakemake --use-conda --cores 32
```

The number of synonymous / non-synonymous sites per gene (`Ls`, `Ln` in
`data/intermediate/raw_genes_filtered_w_theta_hat_*.tsv.gz`) is computed from the prodigal
protein predictions (`combined_reference_*.faa`) by weighting each amino acid by its codon
degeneracy.
