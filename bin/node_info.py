#!/usr/bin/python

# Claatu::node_info -Collect clade info 
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
#  node_info.py  #
#                 #
###################


import dendropy
import os
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("tree_fp", help="file path to the tree that is to be modified")
parser.add_argument("out_fp", help="file path for out files. No extension needed")
parser.add_argument("-p", "--prefix", help = "prefix to file names", default = "claatu")
args = parser.parse_args()

tree_fp = args.tree_fp # should be pickled tree from prep tree
out_fp = args.out_fp
prefix = args.prefix

tree1 = dendropy.Tree.get(path = "{0}".format(tree_fp), schema = "newick")

def GetLevel(tree):
	"return a dict of the level of each node"	
	node_it = tree.preorder_node_iter()
	level_dict = {}
	for node in node_it:
		a = node.level()
		level_dict[str(node.label)] = a
	return level_dict

def	GetDistRoot(tree):
	"return a list of distances from node -> root"
	node_it = tree.preorder_node_iter()
	rootdist_dict = {}
	for node in node_it:
		if node.edge_length == None:
			dist = 0 
			rootdist_dict[str(node.label)] = dist
		else:
			dist = node.distance_from_root()
			rootdist_dict[str(node.label)] = dist
	return rootdist_dict

def Median(list):
	"compute the median value of a list" 
	list = list
	slist = sorted(list)
	list_len = len(slist)
	index = int(list_len/2)
	
	if list_len % 2==0:
		num1 = slist[index]
		num2 = slist[index-1]
		median = float((num1+num2)/2)
	else:
		median = slist[index]
	return median

def GetMedianDist(dist_dict):
	"get the median distance from root to tip derived from all tips"
	vals = []
	for key in dist_dict:
		vals.append(dist_dict[key])
	dict_median = Median(vals)
	return dict_median
		
	
#level_dict = GetLevel(tree1)
#dist_dict = GetDistRoot(tree1)
#dist_median = GetMedianDist(dist_dict)


def DoYourPrint(level_dict, dist_dict, dist_median, out_fp):
	"this function prints out the data to files"
	ldict = "{0}/{1}_{2}".format(out_fp, prefix, "levels.txt")
	ddict = "{0}/{1}_{2}".format(out_fp, prefix, "dist.txt")
	dmedian = "{0}/{1}_{2}".format(out_fp, prefix, "dist_median.txt")
	f1 = open(ldict, 'w+')
	f2 = open(ddict, 'w+')
	f3 = open(dmedian, 'w+')
	#print out node levels
	for item in level_dict:
		print >> f1, "%s\t%s" % (item, level_dict[item])
	f1.close()
	#print out distances
	for dist in dist_dict:
		print >> f2, "%s\t%s" % (dist, dist_dict[dist])
	f2.close()
	#print out median distance from node to root
	print >> f3, dist_median, 	
	f3.close()
	
	return None


level_dict = GetLevel(tree1)
dist_dict = GetDistRoot(tree1)
dist_median = GetMedianDist(dist_dict)	
DoYourPrint(level_dict, dist_dict, dist_median, out_fp)

