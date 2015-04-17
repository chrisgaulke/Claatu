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

'''

This is the beginning of a script that will accomplish the following tasks

1. De novo 16S tree 
	a. Align and construct 16S tree
	b. Construct a names file ( dictionary of rep seqs(keys) and all seqs identical to this seq (vals)  
	c. Import tree
	d. Use tree to make a node x sample file and populate it with data from tree and names file

2. Traversal tool for any ref tree
	a. Names nodes
	b. Makes lookup table of ancestors of each node

3. Read placed data on ref tree
	a. Imports traversal tool names nodes and makes a lookup table. 
	b. Run pplacer, import and parse pplacer outs. 
	c. Collapse placements into nodes and make a samples x nodes count table . 

4. Clade matrix analysis tools 
	a. additional tools for traversing trees

'''
#get working directory. In the future would like to specify 
wd = os.getcwd()

####
#Prep Tree
####

#eventually taken from command line or config file
tree_fp = '/Users/gaulkec/Documents/Coevo_project/Co_evo_test_data/new_test_tree.tre'

#eventually taken from command line or config file
tree_type = "newick"

#might want to save bootstraps for latter
def PrepTree(tree_fp, tree_type):
	#import tree object
	tree1 = dendropy.Tree.get_from_path("{0}".format(tree_fp), schema="{0}".format(tree_type))
	
	#name nodes
	node_it = tree1.preorder_node_iter()
	k = 1
	for i in node_it:
		if i.label == None:
			if hasattr(i, 'taxon') and i.taxon != None:
				i.label = i.taxon.label
			else:
				#continue
				i.label = "root"
		else:
			j = str(k)
			i.label = "{0}{1}".format("node", j) 
			k = k + 1
	return tree1

#write tree so it can be used by other programs
#need to test that this can be used by other programs
tree1.write_to_path(
        'new.tre',
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

####
#External calls
####
''' 
At this point we will need to make an external calls to several programs including pplacer, taxit, FastTree, and infernal.
Should check out 'subprocess.call(['command', '-options'])'.

'''
#Might be best to put the tree pickler, pplacer and .jplace file parser into separate modules.

#after query reads aligned with each other and tree reads then we do so post processing of the data
		
####
#Parse .jplace file
####

# make into a dictionary

#in future this will be read in from previous step automatically
jplace_file = "/Users/gaulkec/Desktop/test.jplace" 

#tested 4/8/15

def PlaceParse(jplace_file):
	parse_dict = {}
	with open (jplace_file, "r") as jfile:
		data=jfile.read().replace('\n', '')
	
	place_elem = data.split('"placements":')
	place_tre = place_elem[0]
	place_else = place_elem[1]
	place_elem = place_else.split('"metadata":')
	place_placements = place_elem[0]
	place_extra = place_elem[1]
	place_elem = place_extra.split('"fields":')
	place_fields = place_elem[1]
	
	parse_dict['tree'] = place_tre
	parse_dict['placements'] = place_placements
	parse_dict['metadata'] = place_elem
	parse_dict['fields'] = place_fields
	
	return parse_dict 

#call	
parse_dict = PlaceParse(jplace_file)

####
#Parse tre string
####
#tree edge lookoup dictionary. Edge ids are the keys node names are the values

#tested 4/8/15
def EdgetoTail(jplace_tree):
	"this function will map all edge labels to their corresponding tail nodes"
	tree_dict = {}
	#parse the tree structure. The edge IDs are contained in a set of {} 
	t1 = jplace_tree.strip('":tree{}; ') #remove from the beginning of the tree 
	t2 = t1.strip(',";') #clean up additional characters on the ends of strings
	l = t2.split("}")#split on the close curly brace which results in exactly 1 edge id /node as the ids always follow immediately after a node
	
	#remove the last element of the list which will be empty
	del l[-1] #should be a ';' as this is the last element of the tree string which follows the last } 
	
	for word in l:
		s1 = word.split("{") 
		edge_id = s1[1]
		s2 = s1[0]
		s3 = s2.split(":") 
		#node name
		name = s3[0]
		name = name.strip(')(,')
		if name == '':
			name = 'root'
			tree_dict[edge_id] = name
		else:
			tree_dict[edge_id] = name 
	return tree_dict

#call	
tree_dict = EdgetoTail(parse_dict['tree'])

####
#Parse Placements
####			

#tested 4/8/15
def Stripper(listx):
	"strips whitespace from each element in a list"
	nlist = []
	for y in listx:
		val = y.replace(' ', '')
		nlist.append(val)
	return nlist 

#tested 4/8/15
def PlaceParser (placement_string):
	"This function will parse all the placements found in the placements slot of the .jplace dictionary"
	#This will initialize a regex pattern for finding items enclosed in []
	pattern = re.compile('\[.*?]')
	place_dict = {}
	place_placements = placement_string
	
	place_placements = place_placements.strip(' [,]')			
	placements = place_placements.split('}')
	placements = filter(None, placements)
	
	for placement in placements: 
		place_count = 0 
		placement = placement.strip('{ "p:,')
	
		placement_stats = placement.split('"nm":')[0].strip(' ,')
		placement_nm = placement.split('"nm":')[1].strip(' ,')
	
		placement_nm = placement_nm[1:-1]
		placement_stats = placement_stats[1:-1]
	
		stats = []
		place = pattern.findall(placement_stats)
		place = filter(None, place)
		place = Stripper(place)
		
		for p in place:
			stats_list = p.strip('[] ').split(',')
			stats.append(stats_list)
			place_count += 1
	
			#make a dictionary of all placements where the key is the read Id
		nm = pattern.findall(placement_nm)
		nm = filter(None, nm)
		for name_mass in nm:
			name = name_mass.split(',')[0].strip('][" ')
			mass = name_mass.split(',')[1].strip('][" ')
			place_dict[name] = [mass, stats] 	
	return place_dict

#call
place_dict = PlaceParser(parse_dict['placements'])


####
#Parse fields
####

#tested 4/8/15
def ParseFields(place_fields):
	"This will parse the fields string from the .jplace file and will return a list of field names in order"
	place_fields = place_fields.strip('{} []')
	stats_fields = place_fields.replace('"', '').split(',')
	fields = []
	for word in stats_fields:
		word = word.strip()
		fields.append(word)	 
	return fields

#call	
fields = ParseFields(parse_dict['fields'])
 	
####
#Parse Metadata
####

#Not sure if we really even need this might push it to a log file in the future. 
''' I am just going to store this as a string until I figure out what to do with it. '''


####
#End .jplace file parser
####

####
#Integration of jplace output and tree data
####


#tested 4/8/15
def MakeTailtoHead(dict):
	"this is a function to make a lookup table of tail nodes to head nodes"
	tail_to_head = {}
	for key, value in dict.items():
		tail_node = tree1.find_node_with_label(value)
		head_node = tail_node.parent_node
		if hasattr(head_node, 'label'):
			tail_to_head[value] = head_node.label
		else:
			continue
	return tail_to_head

#call
tail_to_head_dict = MakeTailtoHead(tree_dict)	
			
	

''' format of place_dict: 

place_dict[name] = [mass, [[stats_0],[stats_1], [stats_n]]]

name: Read/fragment name
mass: mass
stats:["edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length", etc.]

'''

#Need to make a dataset to test this
def GetBest(list):
	"This function gets a best hit for each read and reports a best hits. There are no tie breakers (i.e., the first largest number is counted as the best hit) "
	#need to eventually incorporate a nearest common ancestor function to replace this function
	hit = 0
	hit_name = ''
	for i in list:
		var1 = i[1]
		var2 = i[0]
		if float(var1) > hit:
			hit = float(var1)
			hit_name = var2
		else:
			continue
	return hit_name

#tested 4/8/15
#samples need to be defined by the user or parsed from a file elsewhere but for now...
samples = ['Monkey', 'Human', 'Mouse']  


def MultiplexCountPlaceHits(samples, tree_dict, place_dict):
	"This will do edge -> samps instead of samples -> edges"
	edge_samp = {}
	for edge in tree_dict:
		edge_samp[edge] = {}
	for sample in samples:
		for edge in tree_dict:
			edge_samp[edge][sample] = 0
	for p, s in place_dict.items():
		if len(s) > 2:
			for stat in s[1]:
				best_hit = GetBest(s[1]) #retain only the best hit
				name = p.split('_')[0] #sample ID
				if stat[0] == best_hit and stat[0] in edge_samp:
					edge_samp[stat[0]][name] += 1
				else:
					print 'error {0} not found'.format(stat[0])
		else:
			stat = s[1][0]
			name = p.split('_')[0]
			if stat[0] in edge_samp:
				edge_samp[stat[0]][name] += 1
			else:
				print 'error {0} not found'.format(stat[0])
	return edge_samp

#Call
edge_samps = MultiplexCountPlaceHits(samples, tree_dict, place_dict)


####
#Get Node lookup 
####


#tested 4/9/15
def MultiplexEdgeHeadMapper(edge_samps, tail_to_head_dict, tree_dict):
	"This function will map the edge counts to head nodes"
	head_count_dict = {}
	
	for edge, sdict in edge_samps.items():
		try:
			key = tail_to_head_dict[tree_dict[edge]]
		except KeyError:
			if tree_dict[edge] == 'root':
				key = 'root'
			else:
				print "Hmmm. There seems to be a problem with your keys"
		if key in head_count_dict:
			#sum up all values for all samples
			for entry in sdict:
				my_sum = sdict[entry] + head_count_dict[key][entry]
				head_count_dict[key][entry] = my_sum #test
		else:	
			head_count_dict[key] = sdict
	
	return head_count_dict
	

#Call
head_count_dict = MultiplexEdgeHeadMapper(edge_samps, tail_to_head_dict, tree_dict)

####
#Should start thinking about what I can destroy or pickle at this point
####

#need to collapse this into a lookup table of cumulative

#tested 4/16/15
def CumulativeNodeCounts(head_count_dict, samples):
	"This function will collapse the node counts at each level (i.e., will sum the counts of a n node and it's children and return a list of summed sample counts"
	sum_temp = {}
	
	for node, sdict in head_count_dict.items():
		nodex = tree1.find_node_with_label(node)
		lnodex = nodex.label #label for node
		n_iter = nodex.preorder_internal_node_iter() #make a iterator for the children of each node
		sum_temp[lnodex] = {} #set up a new dictionary
		for sample in samples:
			sum_temp[lnodex][sample] = 0
		for n in n_iter:
			if n.label in head_count_dict:
				for sample in head_count_dict[n.label]:
					sum_temp[lnodex][sample] += head_count_dict[n.label][sample]
			else:
				continue
	return sum_temp

#call 
collapsed_node_counts = CumulativeNodeCounts(head_count_dict, samples)
'''loose ends: 
1) need to figure out how to deal with name files to
	a) get sample names
	b) deal with single sequences that represent multiple seqs
'''

