from Bio import Phylo
import argparse
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd


def DAY(t, t1):
    # Root trees with same root
    root = t.get_terminals()[0].name
    t.root_with_outgroup({"name": root})
    t1.root_with_outgroup({"name": root})
    
    terminals  = [x.name for x in t.get_terminals()]
    terminals1 = [x.name for x in t1.get_terminals()]
    DF         = {}
    
    for i in range(len(terminals)):
        DF[terminals[i]] = i
    
    NODES  = day_nodes(DF, t.clade, [])
    NODES1 = day_nodes(DF, t1.clade, [])

    return NODES, NODES1
    
def day_nodes(DF, c, NODES):
    if c.is_terminal() == False:
        
        l = [DF[x.name] for x in c.get_terminals()]
        l.sort()
        NODES.append(l)
            
        day_nodes(DF, c.clades[0], NODES)
        day_nodes(DF, c.clades[1], NODES)
    
    if len(NODES) == len(c.get_nonterminals()):
        return NODES
    
def RF_dist(tree1, tree2):
    lis = DAY(tree1, tree2)
    score = 0
    for i in lis[0]:
        if i in lis[1]:
            score += 1
    return len(lis[0]) + len(lis[1]) - 2*score


def tree_matrix(path):
	mypath    = path
	trees = [f for f in listdir(mypath) if isfile(join(mypath, f))]
	n_trees = len(trees)
	RF_matrix = np.zeros((n_trees, n_trees))

	for i in range(n_trees):
	    for j in range(n_trees):
	        if i > j and i != j:
	            tree_i = Phylo.read(path+trees[i], "newick")
	            tree_j = Phylo.read(path+trees[j], "newick")
	            RF_matrix[i, j] = RF_dist(tree_i, tree_j)
	            
	RF_matrix = RF_matrix.T + RF_matrix
	
	df = pd.DataFrame(RF_matrix)
	df.columns = trees
	df.index = trees

	return df

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("-t1", help = "Input file 1, format: NEW")
	parser.add_argument("-t2", help = "Input file 2, format: NEW")
	parser.add_argument("-multiple", help = "Find distance between multiple NEW-files in the same directory")
	args   = parser.parse_args()

	if args.multiple != None:
		distances = tree_matrix(args.multiple)
		print(distances)
		return

	tree1  = Phylo.read(args.t1, "newick")
	tree2  = Phylo.read(args.t2, "newick")

	RF     = RF_dist(tree1, tree2)

	print("Distance between trees: ", RF) 


if __name__ == '__main__':
	main()