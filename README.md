ClaaTU
======

A workflow for quantifying clade abundance in phylogenetic trees.

Version: Alpha

Overview
--------

The ClaaTU (**Cla**d**a**l **T**axonomic **U**nits) workflow uses data files produced by third party microbiome analysis software (e.g., QIIME and Mothur) to identify and quantify the abundance of specific clades in a user provided phylogenetic tree. This allows us to examine the abundance of clades across the input tree from tip to root. Currently, ClaaTU is run as a collection of scripts and a workflow is provided below. 

Workflow
--------

1. The first thing that ClaaTU needs is a newick tree without internal node identifiers. In this step ClaaTU decorates internal nodes with node identifiers that will be used in downstream scripts. If internal node identifiers are present in the tree provided Claatu will overwrite them.  

	    python <path_to_ClaaTU/bin/prep_tree.py> <input_newick_tree.tre> 

    The output of this script is a phylogenetic tree named "new_prepped_tree.tre". You will use this tree for the rest of the workflow.    
2. Now that we have the ClaaTU ready tree we can get to work. The next script that we will run will take in a table delimited table of OTU counts. You should have this from you microbiome analysis. Note that ClaaTU will not recognize a biom formatted table, so be sure to convert you biom table to text format. 

	python <path_to_ClaaTU/bin/count_tree.py> <otu_table> <prepped_tree> <out_file_path> 

    The output of this script is a clade counts table with internal node identifiers as columns and sample IDs as rows. These data can be used to examine differential abundance of clades across a case and control study. For some this might be where you stop if all you care about is what clades are in a sample and at what abundance. However, we have add some additional scripts downstream that some might find useful.  

3. We realize that many researchers may be interested in the taxonomic labels associated with each clade, particularly if the clade is significantly different between case and control. Currently, we have attempted to address this need by implementing a strategy of propagation of taxonomic labels that use the taxonomic association of the tips (OTUs) and attempts to assign them to the ancestors of this tip. Importantly, this taxonomic label is only assigned to the ancestor of the tip if 100% of the other tips contained by the parent node also share this label. This script requires a taxonomy file, which may have been an output of upstream microbiome analysis. For example, if you used QIIME to analyze you data this file will be called something like rep_set_tax_assignments.txt. The first two columns of this taxonomy file should be OTU ID followed by a taxonomy string. This file must be tab delimited. Taxonomic label propagation can be performed by executing the following:
 
	python <path_to_ClaaTU/bin/clade_stat.py> <prepped_tree> <tax_file> <out_file_path> -p <file_prefix>

    This generates a tab delimited text file (<file_prefix>_tax2node.txt) with clade ID in the first column and a list of taxonomy strings in the second column.  

4. Next we try to find which taxonomic level (if any) is shared by each tip in a clade. 
	
	```python <path_to_ClaaTU/bin/tax_parser.py> <prefix>_tax2node.txt <outfile_name>```

This will output a final taxonomy dictionary that contains a mapping of OTU (column 1) to tax ID (column 2).

5. Finally, we gather some information about the tree and the nodes in it. 

	```python <path_to_ClaaTU/bin/node_info.py> <prepped_tree> <out_file_path> -p <file_prefix>```

This should produce three files: <file_prefix>_levels.txt which contains information about how nested a node is, <file_prefix>_dist_median.txt which is a single line file containing the median branch length to the root, and <file_prefix>_dist.txt which gives the branch length to the root from each node. 

---

Note
----
Claatu is a work in progress and will be undergoing rapid development during the coming months. Growing pains are to be expected so please let me know if you discover a bug. 
~

