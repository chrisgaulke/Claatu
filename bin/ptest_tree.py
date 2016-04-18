#!/usr/bin/python

# Claatu::count_tree -Construct clade lookup tables 
#Copyright (C) 2016  Christopher A. Gaulke 
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
#  ptest_tree.py  #
#                 #
###################


#for sample random distro

import re
import os
import argparse
import sys
import dendropy
from random import randint
from random import shuffle
from scipy import mean as mean
from numpy import std as sd
from scipy.stats import norm as norm

parser = argparse.ArgumentParser()
parser.add_argument("biom_fp", help="file path to the biom table to be used")
parser.add_argument("tree_fp", help="file path to the tree that is to be modified")
parser.add_argument("out_fp", help="name and file path for out file. Should have .txt or .tab extension")
parser.add_argument("-p", "--perms", help="number of permutations", default = int(1))
parser.add_argument("-m", "--method", help="Method of randomization [samples, labels]", default = "samples")
parser.add_argument("-g", "--groups", help="Mapping of sample ID to group for ptests", default = None)


args = parser.parse_args()
biom_fp = args.biom_fp
out_fp = args.out_fp
perms = int(args.perms)
method = args.method
map = args.groups
tree_fp = args.tree_fp # should be result from prep tree

tree1 = dendropy.Tree.get(path = "{0}".format(tree_fp), schema = "newick")

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

def PermuteBiomLabels(biom_dict):
	"This function will randomly redistribute counts from samples across different otus"
	#create an empty copy of biom_dict
	ndict = {}
	for sample in biom_dict:
		ndict[sample]={}
		for otu in biom_dict[sample]:
			ndict[sample][otu] = 0
	
	for sample in biom_dict:		
		#total = sum(biom_dict[sample].itervalues())
		total=0
		for otu in biom_dict[sample]:
			total += float(biom_dict[sample][otu])
		total = int(total)
		keys = ndict[sample].keys()
		lk = len(keys)-1
		while total > 0:
			mykey = keys[randint(0,lk)]
			ndict[sample][mykey] += 1
			total -= 1  
	return ndict

def PermuteSampleLabels(biom_dict):
	"This function will randomly redistribute count bins from one otu to antoher within a sample"
	#create an empty copy of biom_dict
	ndict = {}
	for sample in biom_dict:
		ndict[sample]={}
		my_vals = biom_dict[sample].values()
		my_vals.reverse()
		shuffle(my_vals)
		for otu in biom_dict[sample]:
			ndict[sample][otu] = my_vals.pop()
	return ndict
			
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
		tip = tip.strip('\'')
		tip_ancestors[tip] = vals
	return tip_ancestors
					
def AncestorCrawl(ancestors, pbiom):
	"This functions collects a list of ancestors of each OTU from a reference tree and assigns the count value of each OTU to the ancestor"
	"""If there are values in the slot already this function will sum the value"""	
	ancestors = ancestors
	#Initialize the dict
	cml_nodes = {}
	for sample in pbiom:
		cml_nodes[sample] = {}
	
	for sample in pbiom:
		for otu in pbiom[sample]:
			if otu in ancestors: #problem begins here
				for ant in ancestors[otu]:
					if ant in cml_nodes[sample]:
						cml_nodes[sample][ant] += float(pbiom[sample][otu])
					else:
						cml_nodes[sample][ant] = 0 
						cml_nodes[sample][ant] += float(pbiom[sample][otu])
	
	return cml_nodes

def CalculateCoreness(cml_nodes):
	"Calculate the coreness of each clade"
	core_dict = {}
	sample_keys = cml_nodes.keys()
	otu_keys = cml_nodes[sample_keys[0]].keys()
	for key in otu_keys:
		core_dict[key] = 0
	for sample in cml_nodes:
		for otu in cml_nodes[sample]:
			if int(cml_nodes[sample][otu]) != 0:
				core_dict[otu] += 1
			else:
				continue
	for key in core_dict.keys():
		core_dict[key] = (float(core_dict[key]) / len(sample_keys))
	return core_dict
			  
def MakeTable(perm_dict, perms, out_fp):
	"this will make a table from combined actual and random coreness data"
	#make a pretty(ish) table 
	my_range = range(1,perms+1)
	my_col_header = ["exp_values"]
	for i in my_range:
		my_col_header.append(i)
	f1 = open(out_fp, 'w+')
	#print header
	print >> f1, "\t".join(str(x) for x in my_col_header),
	print >> f1, "\n",
	#Print out coreness for real data (column1) and randomized data (rest of columns)
	for otu in perm_dict:
		print >> f1, "%s\t" % otu,
		print >> f1, "\t".join("%4.3f" % x for x in perm_dict[otu]),
		print >> f1, "\n",
	
	f1.close()
	return None

def RunIter(biom_fp, perms, tree, out_fp, method): 
	"runs the p-test"
	my_biom = BiomTabParser(biom_fp)
	tip_an_dict = TipAncestorLookup(tree)
	real_cml_node_dict = AncestorCrawl(tip_an_dict, my_biom)
	real_core_dict = CalculateCoreness(real_cml_node_dict)
	#use the real coreness dictionary keys to make a dictionary for permutations
	perm_dict = {}
	nperms = perms
	for otu in real_core_dict:
		perm_dict[otu] = [real_core_dict[otu]]
	#now do the permutations
	if method == "samples":
		while nperms > 0:
			nperms -= 1
			perm_biom = PermuteSampleLabels(my_biom)
			cml_node_dict = AncestorCrawl(tip_an_dict, perm_biom)
			core_dict = CalculateCoreness(cml_node_dict)
			for otu in core_dict:
				perm_dict[otu].append(core_dict[otu])
	else:
		while nperms > 0:
			nperms -= 1
			perm_biom = PermuteBiomLabels(my_biom)
			cml_node_dict = AncestorCrawl(tip_an_dict, perm_biom)
			core_dict = CalculateCoreness(cml_node_dict)
			for otu in core_dict:
				perm_dict[otu].append(core_dict[otu])
	
	MakeTable(perm_dict, perms, out_fp)
	Ztest(perm_dict, out_fp)
	return None 

#need to make a function that will parse a mapping file

def ParseMap(map):
	"This function takes a mapping file and creates a dictionary"
	mdict = {}
	with open(map) as f:
		header = f.readline()
		for line in f:
			line = line.strip()
			elem = line.split('\t')
			#for e in elem:
			#	e = e.rstrip()
			if elem[1] in mdict:
				mdict[elem[1]].append(elem[0])
			else:
				mdict[elem[1]] =[]
				mdict[elem[1]].append(elem[0])				
	return(mdict)


def CalculateGroupCoreness(cml_nodes, map):
	"Calculate the coreness of each clade"
	ldict = {}
	core_dict = {}
	for m in map:
		ldict[m] = len(map[m])
		core_dict[m] ={}
	sample_keys = cml_nodes.keys()
	otu_keys = cml_nodes[sample_keys[0]].keys()
	for group in core_dict:
		for otu in otu_keys:
			core_dict[group][otu] = 0
		for sample in cml_nodes:
			if sample in map[group]:
				for otu in cml_nodes[sample]:
					if int(cml_nodes[sample][otu]) != 0:
						core_dict[group][otu] += 1
					else:
						continue
			else:
				continue
	for group in core_dict:
		for key in core_dict[group].keys():
			core_dict[group][key] = (float(core_dict[group][key]) / ldict[group])
	return core_dict

def RunGroupIter(biom_fp, perms, tree, out_fp, method, map):
	"runs the p-test"
	mmap = ParseMap(map)
	my_biom = BiomTabParser(biom_fp)
	tip_an_dict = TipAncestorLookup(tree)
	real_cml_node_dict = AncestorCrawl(tip_an_dict, my_biom)
	real_core_dict = CalculateGroupCoreness(real_cml_node_dict, mmap)
	#use the real coreness dictionary keys to make a dictionary for permutations
	perm_dict = {}
	nperms = perms
	
	for group in real_core_dict:
		for clade in real_core_dict[group]:
			perm_dict['_'.join([clade,group])] = [real_core_dict[group][clade]]
	
	#now do the permutations
	if method == "samples":
		while nperms > 0:
			nperms -= 1
			perm_biom = PermuteSampleLabels(my_biom)
			cml_node_dict = AncestorCrawl(tip_an_dict, perm_biom)
			core_dict = CalculateGroupCoreness(cml_node_dict, mmap)
			for group in core_dict:
				for clade in core_dict[group]:
					perm_dict['_'.join([clade,group])].append(core_dict[group][clade])
	else:
		while nperms > 0:
			nperms -= 1
			perm_biom = PermuteBiomLabels(my_biom)
			cml_node_dict = AncestorCrawl(tip_an_dict, perm_biom)
			core_dict = CalculateGroupCoreness(cml_node_dict, mmap)
			for group in core_dict:
				for clade in core_dict[group]:
					perm_dict['_'.join([clade,group])].append(core_dict[group][clade])
	MakeTable(perm_dict, perms, out_fp)
	Ztest(perm_dict, out_fp)
	return None	 

def Ztest(perm_dict, out_fp):
	"this calculate a zscore and zstat from perm_dict values"
	#Print out ztest results
	#Columns are 1) exp coreness 2) Permutation mean 3) Perm sd 4) clade zscore 5) pval
	out_file = out_fp + "_stats.txt" 
	f1 = open(out_file, 'w+')
	for otu in perm_dict:
		exp = perm_dict[otu][0] # observed coreness
		m = mean(perm_dict[otu][1:]) # mean perms 
		s = sd(perm_dict[otu][1:]) # standard deviation perms
		if s == 0: # might need to develop a different solution for this. 
			z = 1
			p = norm.sf(z) #upper tail of cumulative probability distribution
			print >> f1, "%s\t" % otu,
			print >> f1, "\t".join("%.2E" % x for x in [exp,m,s,z,p]), #join stats and print
			print >> f1, "\n",	
		else: 
			z = (float(exp) - float(m)) / float(s) #zscore
			p = norm.sf(z) #upper tail of cumulative probability distribution
			print >> f1, "%s\t" % otu,
			#print >> f1, "\t".join("%4.3f" % x for x in [exp,m,s,z,p]), #join stats and print
			print >> f1, "\t".join("%.2E" % x for x in [exp,m,s,z,p]), #join stats and print
			print >> f1, "\n",	
	f1.close()
	return None


if map == None:
	RunIter(biom_fp, perms, tree1, out_fp, method)
else:
	RunGroupIter(biom_fp, perms, tree1, out_fp, method, map)
	

