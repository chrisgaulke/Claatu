#!/usr/bin/python

# Claatu::claatu_calcs - Tests for nestedness 
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

#####################
#                   #
#  claatu_calcs.py  #
#                   #
#####################


import dendropy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("tree_fp", help="file path to the tree that is to be modified")
parser.add_argument("node_fp", help="file path to the a list of nodes (one per line) that are of interest")
args = parser.parse_args()
tree_fp = args.tree_fp
node_fp = args.node_fp

tree1 = dendropy.Tree.get(path = "{0}".format(tree_fp), schema = "newick")

def AncestorLookup(tree):
	"This function makes a dictionary of the ancestors of each node"
	node_it = tree.preorder_node_iter()	
	tip_ancestors = {}
	#make an iterator for each node and append each ancestor to a list(vals)
	for node in node_it:
		ancest_it = node.ancestor_iter(inclusive=False) #get iter for all ancestors	
		vals = []
		for ancestor in ancest_it: 		 
			vals.append(str(ancestor.label))
		tip_ancestors[str(node.label)] = vals
	return tip_ancestors

def MakeNodeList(node_file):
	"This func takes a file with one node per line and makes into a list for IsNested" 
	with open(node_file) as f1:
		nl = []
		for line in f1:
			line.strip()
			nl.append(line)
	return nl

#This consumes lots of ram and might not be the right way to go about this. 
#So I will only do it for the nodes that are interesting
 
def IsNested(node_list, an_dict):
	"Test if one node is nested in another"
	nest_dict ={}
	#node_list = node_list
	for nodeA in node_list:
		nodeA = nodeA.strip() 
		nest_dict[nodeA] = {}
		for nodeB in node_list:
			nodeB = nodeB.strip()
			if nodeB in an_dict[nodeA]:
				nest_dict[nodeA][nodeB] = "True"
			else:
				nest_dict[nodeA][nodeB] = "False"
	return nest_dict
	

an_dict = AncestorLookup(tree1)
nl = MakeNodeList(node_fp)
nest_dict = IsNested(nl, an_dict)
print nest_dict