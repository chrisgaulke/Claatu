#!/usr/bin/python

# Claatu::count_tree -Construct clade lookup tables 
#Copyright (C) 2015  Christopher A. Gaulke 
#author contact: gaulkec@science.oregonstate.edu
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#    
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#    
#You should have received a copy of the GNU General Public License
#along with this program (see LICENSE.txt).  If not, see 
#<http://www.gnu.org/licenses/>

####################
#   ___________    #
#  |           |   #
#  | |==(*)==| |   #
#   |         |    #
#    |_______|     #
#                  #
####################  

###################
#                 #
#  count_tree.py  #
#                 #
###################


import dendropy
import re
import os
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("biom_fp", help="file path to the biom table to be used")
parser.add_argument("tree_fp", help="file path to the tree that is to be modified")
parser.add_argument("out_fp", help="file path for out file. Should have .txt or .tab extension")
args = parser.parse_args()

biom_fp = args.biom_fp
tree_fp = args.tree_fp # should be result from prep tree
out_fp = args.out_fp

tree1 = dendropy.Tree.get(path = "{0}".format(tree_fp), schema = "newick")
#tree1 = pickle.load(open("/Users/gaulkec/Desktop/Claatu_pickles/prepped_tre.pkl", "rb"))


def BiomTabParser(biom):
	"This function parses a biom file (closed ref OTUS) in tab delim format and creates a dictionary of samples by OTU by count"
	"""Columns of the table should be samples and rows should be OTUs
	should ignore lines until matching a line #OTU"""
	#set pattern	
	my_pattern = re.compile('\#OTU')
	#read in biom table line by line
	biom_dict = {} 
	my_header = ""
	with open(biom) as f:
		for line in f:
			if re.match('#', line):
				if re.match(my_pattern, line):
					my_header = line
					my_header = my_header.strip('#')
					my_header = my_header.rstrip()
					my_header = my_header.split('\t')
					for i in range(1,len(my_header)):
						biom_dict[my_header[i]] = {}
			else:
				line = line.rstrip()
				temp = line.split("\t")
				for i in range(1,len(temp)):
					biom_dict[my_header[i]][temp[0]] = temp[i]
				
	return biom_dict

def TipAncestorLookup(tree):
	"This function makes a dictionary of the ancestors of each node"
	node_it = tree.leaf_node_iter()
	tip_ancestors = {}
	#make an iterator for each node and append each ancestor to a list(vals)
	for node in node_it:
		ancest_it = node.ancestor_iter(inclusive=False) #get iter for all ancestors	
		vals = []
		for ancestor in ancest_it:
			vals.append(str(ancestor.label))
		tip = str(node.taxon)
		#print tip
		tip = tip.strip('\'')
		#print tip
		#tip_ancestors[str(node.taxon)] = vals
		tip_ancestors[tip] = vals
	return tip_ancestors
					

def AncestorCrawl(ancestors, pbiom):
	"This functions collects a list of ancestors of each OTU from a reference tree and assigns the count value of each OTU to the ancestor"
	"""If there are values in the slot already this function will sum the value"""	
	ancestors = ancestors
	#print ancestors
	#Initialize the dict
	cml_nodes = {}
	for sample in pbiom:
		cml_nodes[sample] = {}
	
	for sample in pbiom:
		for otu in pbiom[sample]:
			#print otu
			if otu in ancestors: #problem begins here
				for ant in ancestors[otu]:
					#print ant
					if ant in cml_nodes[sample]:
						cml_nodes[sample][ant] += float(pbiom[sample][otu])
					else:
						cml_nodes[sample][ant] = 0 
						cml_nodes[sample][ant] += float(pbiom[sample][otu])
	
	return cml_nodes
			  
def MakeTable(cml_node_dict, out_fp):
	"this will make a table from cml_node_dict"
	#make a pretty(ish) table 
	my_col_header = ["Sample"]
	my_row_header =[]	
	#sample list
	for sample in cml_node_dict:
		my_row_header.append(str(sample))
	#OTU list
	for otu in cml_node_dict[my_row_header[0]]:
		my_col_header.append(str(otu))
	#make file
	f1 = open(out_fp, 'w+')
	#print header
	print >> f1, "\t".join(my_col_header),
	print >> f1, "\n",
	#Print sample by OTUs
	for sample in my_row_header:
		print >> f1,"%s\t" % (sample),
		for i in range(1, len(my_col_header)):
			print >> f1, "%s\t" % (cml_node_dict[sample][my_col_header[i]]),
		print >> f1, "\n",
	f1.close()
	return None


my_biom = BiomTabParser(biom_fp)
tip_an_dict = TipAncestorLookup(tree1)
#print tip_an_dict
cml_node_dict = AncestorCrawl(tip_an_dict, my_biom)	
#print cml_node_dict
MakeTable(cml_node_dict, out_fp)	

