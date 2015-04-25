#!/usr/bin/python

#This script is the beginning of the Claatu software package 


#################
#			    #
#  travtree.py  #
#               #
#################


import dendropy
import pickle
import re
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("tree_fp", help="file path to the tree that is to be modified")
args = parser.parse_args()

#get working directory. In the future would like to specify 
wd = os.getcwd()

####
#Prep Tree
####

#eventually taken from command line or config file
#if you use new.tre you must reset node1 to root bc it is already named. 
#tree_fp = '/Users/gaulkec/Desktop/Claatu_test_files/new.tre'
tree_fp = args.tree_fp

#eventually taken from command line or config file
tree_type = "newick"

#might want to save bootstraps for latter
#this labels tips as tip taxon (i.e., OTU or species name)

def PrepTree(tree_fp, tree_type):
	#import tree object
	tree1 = dendropy.Tree.get_from_path("{0}".format(tree_fp), schema="{0}".format(tree_type))
	
	#name nodes
	node_it = tree1.preorder_node_iter()
	k = 1
	for i in node_it:
		if i.label == None:
			if hasattr(i, 'taxon') and i.taxon != None: # (i.e., a tip)
				i.label = i.taxon.label
			else:
				#continue
				i.label = "root"
		else:
			j = str(k)
			i.label = "{0}{1}".format("node", j) 
			k = k + 1
	return tree1



####
#Make node ancestor lookup table
####

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



def PickleMeTimbers(prepped_tree_obj, ancestor_lookup_table):
	"This function will make a directory for pickled files and pickle the tree and lookup for latter use"
	
	#dir name Claatu_pickles
	pd = '{0}{1}'.format(wd, '/Claatu_pickles')
	try:
		os.stat(pd)
	except:
		os.makedirs('{0}{1}'.format(wd, '/Claatu_pickles'))
	#pickle the tree
	pickled_tree = open('{0}/prepped_tre.pkl'.format(pd), 'wb')
	pickle.dump(prepped_tree_obj, pickled_tree)
	pickled_tree.close()
	
	#pickle lookup table
	pickled_lookup = open('{0}/ancestor_lookup.pkl'.format(pd), 'wb')
	pickle.dump(ancestor_lookup_table, pickled_lookup)
	pickled_lookup.close()
	return None

tree1 = PrepTree(tree_fp, tree_type)  	

ancestor_lookup_dict = AncestorLookup(tree1)

PickleMeTimbers(tree1, ancestor_lookup_dict)

#write tree so it can be used by other programs
#need to test that this can be used by other programs
tree1.write_to_path(
        'new_prepped_tree.tre',
        'newick',
        taxon_set=None,
        suppress_leaf_taxon_labels=False,
        suppress_leaf_node_labels=True,
        suppress_internal_taxon_labels=False,
        suppress_internal_node_labels=False,
        suppress_rooting=False,
        suppress_edge_lengths=False,
        unquoted_underscores=False,
        preserve_spaces=False,
        store_tree_weights=False,
        suppress_annotations=True,
        annotations_as_nhx=False,
        suppress_item_comments=True,
        node_label_element_separator=' ',
        node_label_compose_func=None)