# Study of human fibulin protein of 8 species 

## Uses of the advance_usage.py script: advanced 

Creation of fasta files containing the homologous and paralogous sequences of the two proteique family fibulin and LTBP (Latent Transforming growth factor binding proteins, the outgroup). As well as the gene tree for these two family. 

Script : advanced_usage.sh

Required file : 
- fuse_phylotree.sif : image singularity 
- human_fibulins_ids.txt : refseq identifiers for all human fibulin isoforms 
- human_outgroup_ids.txt : refseq identifiers for all human LTBP isoforms 

File required for anlayse : 
- fibulin_60.fasta 
- no_root_tree.tree 

## Creating the gene tree independently : make_tree_clean

After checking the protein sequences selected with the advance_usage.py script, one sequence was deleted as too small. The gene tree was recalculated using the fasta file fibulin_outgroup.fasta without this sequence. 

Script : make_tree.sh

Required file : 
- 8_species.tree
- fibuline_outgroup_clean.fasta

File required for anlayse : 
- gene_phylo_dir_fibuline_outgroup_clean/clean.tree


## Fuse-phylotree launch on fibulins : run_singularity_fibulin

Script : lanceur_singularity.sh 

Input : 
- fibulin_59.fasta
- fibulin_tree_root.tree
- ppi_shared.csv 

The Itol tree corresponding these results is available [Here](https://itol.embl.de/tree/13125423118220571734007695).
