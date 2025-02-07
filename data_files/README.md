Here's all the files you need to repeat my analyses

Really you only need raw_genes_filtered_w_theta.tsv (description below, under Other files)

All files use the same naming convention more or less so here's a breakdown of of the short hand:

IS -> InStrain (the SNP & analysis pipeline, https://instrain.readthedocs.io/en/latest/)<br>
M[1/5] -> the sample sites M1 or M5<br>
M[1/5]_[3/25] -> the sample site plus MAG origin (MAGs were generated at 3cm below sea floor (cmbsf) and 25cmbsf)<br>
M[1/5]_[3/25]_[1-25] -> Sample site, MAG origin, and read origin depth (in cmbsf)
raw -> not normallized<br>
proc -> filtered (according to IS standards) and normallized (w/ https://github.com/Winshipe/Normalise-Paired-Regions)<br>
stb -> scaffold to bin, associates contig names w/ the MAGs<br>


Other files:<br>
raw_genes_filtered_w_theta.tsv -> This file contains the instrain output for each gene plus theta-hat-n and theta-hat-s and a slightly different Pn & Ps <br>
combined_annotations.txt -> Eggnog (https://github.com/eggnogdb/eggnog-mapper) output<br>
clustered_aaai50_cluster.tsv -> genes clustered w/ mmseqs2 easy-cluster at 50% amino acid identity (if you've got something that really works to cluster these things in a meaningful way please message me on linkedin!)<br>
file2family.txt -> limited taxonomic information from gtdb<br>


