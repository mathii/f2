# Runs on (slighty processed) macs output files to extract the real f2 haplotypes. 

from __future__ import division
import numpy as np
import gzip, sys, getopt
from Bio import Phylo
import pdb

EPS1=1e-4                                  # small number to compare times 
EPS2=1e10                                     # assume that this implies a tree change

##########################################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "tree_file":"", "tree_len_file":"", "out_file":"", "tree_pos_file":"", "resolution":1, "verbose":False  }

    try:
        opts, args = getopt.getopt(sys.argv[1:], "t:p:l:o:r:v:", ["gt", "p", "len", "out", "res", "verbose"])
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        if o in ["-t","--tree"]:      options["tree_file"] = a
        elif o in ["-l","--len"]:     options["tree_len_file"] = a
        elif o in ["-o","--out"]:     options["out_file"] = a
        elif o in ["-p","--pos"]:     options["tree_pos_file"] = a
        elif o in ["-r","--res"]:     options["resolution"] = int(a)
        elif o in ["-v","--ver"]:     options["verbose"] = bool(a)

    print "found options:"
    print options

    return options

##########################################################################################################

def main(options):
    trees_file=gzip.open(options["tree_file"], "r")
    trees=Phylo.parse(trees_file, "newick")

    # load tree positions. 
    positions=gzip.open(options["tree_len_file"], "r")
    pos=positions.read().splitlines()
    pos=np.array([0]+[int(x) for x in pos])
    pos=np.cumsum(pos)

    out=gzip.open(options["out_file"], "w")
    out.write("ID1\tID2\tStart\tEnd\tAge\n")
    
    i=0
    chunks=0
    N=None
    chunk_list=None
    last_pos=-100000
    
    for tree in trees:

        if pos[i]-last_pos<options["resolution"]:
            i+=1
            continue
        else:
            last_pos=pos[i]

        if options["verbose"]:
            print "\rTree " + str(i) +", pos "+ str(pos[i]) + ", " + str(chunks)+" chunks",
            sys.stdout.flush()
        
        if None==N:
            N=tree.count_terminals()
            last_change=np.zeros((N,N), np.int)   # The last position that the TMRCA of i and j changed
            is_f2=np.zeros((N,N), np.bool)         # True if i and j were a cherry in for their last TMRCA 
            f2_age=np.zeros((N,N), np.float64)
            tmrca=np.zeros((N,N), np.float64)         # TMRCA for i and j
            this_tmrca, this_subclade=find_tmrca_and_subclade(tree)
            last_tree=tree
            

        # First find all the nearest neighbour nodes on this tree
        terminal_nodes=tree.get_terminals()
        internal_nodes=tree.get_nonterminals()

        last_tmrca, last_subclade=this_tmrca, this_subclade
        this_tmrca, this_subclade=find_tmrca_and_subclade(tree)

        for node_1 in terminal_nodes:
            label_1=int(node_1.confidence)
            for node_2 in terminal_nodes:
                label_2=int(node_2.confidence)
                if label_1>=label_2:
                    continue

                # this_tmrca=tree.distance(node_1, node_2)/2
                tmrca_change=np.abs(this_tmrca[label_1,label_2]-last_tmrca[label_1, label_2])

                if tmrca_change>EPS1:
                    if is_f2[label_1, label_2]:
                        write_line([label_1, label_2, last_change[label_1,label_2], pos[i], f2_age[label_1,label_2]], out, options)
                        chunks+=1           
  
                    last_change[label_1,label_2]=pos[i]
                    is_f2[label_1,label_2]=False
                    f2_age[label_1,label_2]=-1
                    
        # now, update cherries. 
        for node in internal_nodes:
            if node.is_preterminal(): 
                labels=[int(x.confidence) for x in node.get_terminals()] 
                labels.sort()
                is_f2[labels[0],labels[1]]=True
                f2_age[labels[0],labels[1]]=this_tmrca[labels[0],labels[1]]
                
        last_tree=tree
        i+=1

    ## after the last position, output everything left
    for label_1 in range(N):
        for label_2 in range(label_1, N):
            if is_f2[label_1,label_2]:
                write_line([label_1, label_2, last_change[label_1,label_2], pos[-1], f2_age[label_1,label_2]], out, options)
                chunks+=1


##########################################################################################################

def clade_terminals_changed(tree, last_tree, label_1, label_2):
    """
    Consider the set of terminal nodes of the minimal subclade containing both 
    label1 and label2. Are they the same for tree and last tree? - i.e. did we just 
    graft the node1,node2 subclade somwehere else?
    """
    last_node_1=None
    last_node_2=None
    node_1=None
    node_2=None

    for node in last_tree.get_terminals():
        if int(node.confidence)==label_1:
            last_node_1=node
        elif int(node.confidence)==label_2:
            last_node_2=node
            
    for node in tree.get_terminals():
        if int(node.confidence)==label_1:
            node_1=node
        elif int(node.confidence)==label_2:
            node_2=node

    last_subclade_nodes=last_tree.common_ancestor(last_node_1,last_node_2).get_terminals()
    last_subclade_labels=set([int(x.confidence) for x in last_subclade_nodes])
    this_subclade_nodes=tree.common_ancestor(node_1,node_2).get_terminals()
    this_subclade_labels=set([int(x.confidence) for x in this_subclade_nodes])
    return last_subclade_labels != this_subclade_labels
    
##########################################################################################################

def find_tmrca_and_subclade(tree):
    """
    This finds both the tmrca of each pair of nodes, and the set of node labels which 
    are in the subclade of each pair. i.e. all the terminal nodes of the mrca of the 
    pair of nodes
    """
    N=tree.count_terminals()

    # Initialise tmrca
    tmrca=np.zeros((N,N), np.float64)
    total_depth=max(tree.depths().values())
    tmrca.fill(total_depth)

    # initialise mrca
    subclade_labels={}
    for i in range(N):
        for j in range(i,N):
            subclade_labels[(i,j)]=set(range(N))    

    for node in tree.get_nonterminals():
        leaves=node.get_terminals()
        labels=set([int(x.confidence) for x in leaves] )
        if len(leaves)==N:              # first(?) node is everything. 
            continue
        else:
            for i in labels:
                for j in labels:
                    tmrca[i,j]-=node.branch_length
                    if i<j and labels < subclade_labels[(i,j)]:
                        subclade_labels[(i,j)]=labels
                    
    return tmrca, subclade_labels
           
##########################################################################################################

def write_line(line, out, options):
    if line[3]-line[2]>0:          # if the chunk is >0 length
        output=[int(line[0]/2)+1, int(line[1]/2)+1, int(line[2]), int(line[3]), line[4]]
        out.write("\t".join([str(x) for x in output])+"\n") 

##########################################################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

                       
        
