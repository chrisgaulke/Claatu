#!/usr/bin/python

# Claatu::tax_parser -Gather taxonomic info for nodes
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
#   tax_parser.py   #
#                   #
#####################


#import dendropy
import argparse

#might want ot role this into other scripts

parser = argparse.ArgumentParser()
parser.add_argument("tax_fp", help="file path to the tax dict")
parser.add_argument("out_fp", help="out file path")

args = parser.parse_args()
tax_fp = args.tax_fp
out_fp = args.out_fp
#def ParseNode(node):


def MakeDict(tax_fp):
	" This function is intended to count up all the instances in a tax_dictionary"
	with open(tax_fp) as f1:
		dict = {}
		for line in f1:
			line = line.strip()
			lines = line.split('\t')
			node = lines[0]
			dict[node] = {'kingdom':{}, 'phylum':{}, 'class':{}, 'order':{}, 'family':{}, 'genus':{}, 'species':{}}
			tax = lines[1]
			tax = tax.strip().split(',')
			#print tax
		
			for otu in tax:
				if otu == '':
					continue
				otu = otu.strip()
				divs = otu.split(';')
				bkingdom = divs[0].strip()
				bphylum  = divs[1].strip()
				bclass   = divs[2].strip()
				border   = divs[3].strip()
				bfamily  = divs[4].strip()
				bgenus   = divs[5].strip()
				bspecies = divs[6].strip()
				
				if bkingdom != '' and bkingdom in dict[node]['kingdom']:
					dict[node]['kingdom'][bkingdom] += 1
				elif bkingdom == '' or bkingdom == None:
					continue
				else:
					dict[node]['kingdom'][bkingdom] = 1 
				if bphylum != '' and bphylum in dict[node]['phylum']:
					dict[node]['phylum'][bphylum] += 1
				elif bphylum == '' or bphylum == None:
					continue
				else:
					dict[node]['phylum'][bphylum] = 1 
				if bclass != '' and bclass in dict[node]['class']:
					dict[node]['class'][bclass] += 1
				elif bclass == '' or bclass == None:
					continue
				else:
					dict[node]['class'][bclass] = 1 
				if border !='' and border in dict[node]['order']:
					dict[node]['order'][border] += 1
				elif border == '' or border == None:
					continue
				else:
					dict[node]['order'][border] = 1
				if bfamily !='' and bfamily in dict[node]['family']:
					dict[node]['family'][bfamily] += 1
				elif bfamily == '' or bfamily == None:
					continue
				else:
					dict[node]['family'][bfamily] = 1
				if bgenus !='' and bgenus in dict[node]['genus']:
					dict[node]['genus'][bgenus] += 1
				elif bgenus == '' or bgenus == None:
					continue
				else:
					dict[node]['genus'][bgenus] = 1
				if bspecies !='' and bspecies in dict[node]['species']:
					dict[node]['species'][bspecies] += 1
				elif bspecies == '' or bspecies == None:
					continue
				else:
					dict[node]['species'][bspecies] = 1
				
	return dict


def GetNodeTax(tax_dict):
	"This alogorithm attempts to find the the lowest tax identification (i.e., species, genus, etc.) for each node"
	f_dict ={}
	tax_labels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
	for node in tax_dict:
		if len(tax_dict[node]['kingdom']) > 1:
			f_dict[node] = "mixed_bacteria_archaea"
			continue 
		for label in tax_dict[node]:
			if len(tax_dict[node][tax_labels[0]]) == len(tax_dict[node][tax_labels[6]]) and tax_dict[node][tax_labels[6]].keys()[0] != 's__' :
				for species in tax_dict[node][tax_labels[6]]:
					for genus in tax_dict[node][tax_labels[5]]:
						f_dict[node] = (genus + ';' + species)
			elif len(tax_dict[node][tax_labels[0]]) == len(tax_dict[node][tax_labels[5]]) and tax_dict[node][tax_labels[5]].keys()[0] != 'g__':
				for key in tax_dict[node][tax_labels[5]]:
					f_dict[node] = key
			elif len(tax_dict[node][tax_labels[0]]) == len(tax_dict[node][tax_labels[4]]) and tax_dict[node][tax_labels[4]].keys()[0] != 'f__':
				for key in tax_dict[node][tax_labels[4]]:
					f_dict[node] = key
			elif len(tax_dict[node][tax_labels[0]]) == len(tax_dict[node][tax_labels[3]]) and tax_dict[node][tax_labels[3]].keys()[0] != 'o__':
				for key in tax_dict[node][tax_labels[3]]:
					f_dict[node] = key		
			elif len(tax_dict[node][tax_labels[0]]) == len(tax_dict[node][tax_labels[2]]) and tax_dict[node][tax_labels[2]].keys()[0] != 'c__':
				for key in tax_dict[node][tax_labels[2]]:
					f_dict[node] = key			
			elif len(tax_dict[node][tax_labels[0]]) == len(tax_dict[node][tax_labels[1]]) and tax_dict[node][tax_labels[1]].keys()[0] != 'p__':
				for key in tax_dict[node][tax_labels[1]]:
					f_dict[node] = key
			else:
				for key in tax_dict[node][tax_labels[0]]:
					f_dict[node] = key
	return f_dict
			
def PrintYourStuff(ftax_dict):
	"Prints stuff"
	with open(out_fp, 'w+') as f1: 
		for node in ftax_dict:
			print >> f1, '%s\t%s' % ( node, ftax_dict[node])	 
	return None
	

tax_dict = MakeDict(tax_fp)					
final_tax_dict = GetNodeTax(tax_dict)
PrintYourStuff(final_tax_dict)	
					