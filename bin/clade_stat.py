#!/usr/bin/python

#This script is the part of the Claatu software package 
#These functions might be moved in clade_stat in the future


###################
#                 #
#  clade_stat.py  #
#                 #
###################


import dendropy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("tree_fp", help="file path to the tree that is to be modified")
parser.add_argument("tax_fp", help="file path to the file that contains otus and their tax strings")
parser.add_argument("out_fp", help="file path for out file. Should have .txt or .tab extension")
args = parser.parse_args()

tree_fp = args.tree_fp # should be result from prep tree
tax_fp = args.tax_fp
out_fp = args.out_fp

tree1 = dendropy.Tree.get(path = "{0}".format(tree_fp), schema = "newick")


def GetCladeSizes(tree):
	"This function finds the size of each clade for clade size normalization"
	node_it = tree.preorder_internal_node_iter()
	size_dict = {}
	for node in node_it:
		count = 0
		n_it = node.postorder_iter()
		for n in n_it:
			count += 1
		size_dict[node.label] = count
	return size_dict

def AssignOTULabels2Nodes(tree):
	"This function will assign tax labels to each node" 
	"Gets tax info from each tip and then assigns to ancestors" 
	node_it = tree.preorder_internal_node_iter()
	tips_dict = {}
	for node in node_it:
		tips = node.leaf_nodes()
		vals =[]
		for tip in tips:
			tip = str(tip.taxon)
			tip = tip.strip('\'')
			vals.append(tip)
		tips_dict[node.label] = vals
	return tips_dict

def BuildTaxDict(tax_fp):
	"makes a dictionary of greengenes tax file"
	tax_dict = {}
	with open(tax_fp) as file:
		for line in file:
			entry = line.split('\t')
			otu = entry[0]
			tax_string = entry[1]
			tax_string = tax_string.strip()
			tax_dict[otu] = tax_string
	return tax_dict

def MapTax2Nodes(tax, node_map):
	"This function maps tax strings on the otu IDs and returns a new dictionary"
	otu2tax_dict = {}
	for node in node_map:
		vals =[]
		for otu in node_map[node]:
			vals.append(tax[otu])
		otu2tax_dict[node] = vals
	return otu2tax_dict
	
def WriteFiles(clade_size, ftax, out_fp):
	"this will make a table from cml_node_dict"
	ftax_lab = "{0}_{1}".format(out_fp, "nodes2tax.txt")
	clade_size_lab = "{0}_{1}".format(out_fp, "clade_size.txt")

	f1 = open(ftax_lab, 'w+')
	for node in ftax:
		print >> f1, "%s\t" % (node),
		for tax in ftax[node]:
			print >> f1, "%s," %(tax),
		print >> f1, "\n",
	f1.close()
	f2 = open(clade_size_lab, 'w+')
	for node in clade_size:
		print >> f2, "%s\t%s" % (node, clade_size[node])
	f2.close()
	return None


clade_sizes = GetCladeSizes(tree1)	
#print clade_sizes
otu2tax = AssignOTULabels2Nodes(tree1)	
#print otu2tax		
tax_dict = BuildTaxDict(tax_fp)
#print tax_dict
big_dict = MapTax2Nodes(tax_dict, otu2tax)
#print big_dict
WriteFiles(clade_sizes, big_dict, out_fp)
