#!/usr/bin/python

# Claatu::cut_trees -make subtrees and node lists 
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
#   cut_trees.py    #
#                   #
#####################


import argparse
import dendropy

parser = argparse.ArgumentParser()
parser.add_argument("tree_fp", help="file path to the tree that is to be modified")
parser.add_argument("node_list_fp", help="file path to the node list used for making subtrees")
parser.add_argument("out_fp", help="output tree path")
args = parser.parse_args()

tree_fp = args.tree_fp # should be result from prep tree
node_list_fp = args.node_list_fp
out_fp = args.out_fp

def MakeSubtree(tree_fp, node, out_fp):
	"This will take a tree and a node object as an arg and split a tree making the node the need seed node"
	tree1 = dendropy.Tree.get(path = "{0}".format(tree_fp), schema = "newick")
	node1 = tree1.find_node_with_label(node)
	tree1.seed_node = node1
	tree1.seed_node.parent_node = None
	tree1.update_bipartitions()
	tree1.write(
		path = "%s%s" % (out_fp, node),
		schema= "newick"
		)
	node_it = tree1.preorder_node_iter()
	f2 = "%s%s%s" % (out_fp, node, "nodes.txt")
	with open(f2, 'w+') as file2:
		for node in node_it:
			if node.is_internal():
				print >> file2, "%s" % (node.label)
			elif node.is_leaf():
				tip = str(node.taxon).strip('\'')
				print >> file2, "%s" % (tip)
			else:
				print >> file2, "%s" % ("NA")
				
	return None

def StartCutting(node_list_fp):
	'This function takes a list of nodes and makes subtrees for each of them'
	node_list = []
	with open(node_list_fp) as f1:
		for line in f1:
			line = line.strip()
			node_list.append(line)

	for node in node_list:
		node = node.strip()
		MakeSubtree(tree_fp, node, out_fp)
	return None


StartCutting(node_list_fp)

