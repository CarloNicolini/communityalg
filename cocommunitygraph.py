#!/usr/bin/python

import sys
import numpy as np
import networkx as nx
from scipy.sparse import csr_matrix as np_to_scipy_matrix
import matplotlib as mpl
import matplotlib.pyplot as plt
import random
import find_intersections as fi

def writedot(A,filename):
    graphstr = "graph G{\n"
    graphoptions = "node [style = rounded];\n\
    node [style = \"\"];\n"

    #print [ 'B'+str(j) for j in range(1,A.shape[1])]
    commsastr = "{rank=same " + str([ 'A'+str(j) for j in range(1,A.shape[0])]).replace('[','').replace(']','').replace('\'','')+"}\n"
    commsbstr = "{rank=same " + str([ 'B'+str(j) for j in range(1,A.shape[1])]).replace('[','').replace(']','').replace('\'','')+"}\n"
    
    f=open(filename + '.dot','w')
    f.write(graphstr)
    f.write(graphoptions)
    f.write(commsastr)
    f.write(commsbstr)

    for i in range(0,A.shape[0]):
        for j in range(0,A.shape[1]):
            if A[i,j]!=0:
                s = ('A%d--B%d [penwidth=%1.1f];\n') % (i,j,np.log(A[i,j]+1))
                f.write(s)
    f.write('}\n')

def writepajek(A,filename):
    SP = np_to_scipy_matrix(A)
    edges = []
    for i in range(0,A.shape[0]):
        for j in range(0,A.shape[1]):
            if A[i,j]!=0:
                edges.append(['A'+str(i),'B'+str(j)])

    G = nx.bipartite.from_biadjacency_matrix(SP,create_using=None)
    nx.write_pajek(G, filename + '.net')

def drawbipartite(A):
    # Create the nx bipartite graph
    SP = np_to_scipy_matrix(A)

    B = nx.bipartite.from_biadjacency_matrix(SP,create_using=None)
    m = nx.number_of_edges(B)
    n = nx.number_of_nodes(B)
    B  = nx.convert_node_labels_to_integers(B)
    #B = nx.relabel_nodes(B,lambda x: x+1)
    # Compute the nodes positions
    nodesleft,nodesright = nx.bipartite.sets(B)
    n1,n2 = len(nodesleft),len(nodesright)
    nLeft = len(nodesleft)
    nRight = len(nodesright)

    nodepos = dict()
    nodepos.update( [n,(-1,max(nLeft,nRight)-i) ] for i,n in enumerate(nodesleft))
    nodepos.update( [n,(1,max(nLeft,nRight)-i) ] for i,n in enumerate(nodesright))
    
    # Do the node position permutation
    lines = []
    for e in B.edges():
        x0 = nodepos[e[0]][0]
        y0 = nodepos[e[0]][1]
        x1 = nodepos[e[1]][0]
        y1 = nodepos[e[1]][1]
        lines.append( np.array([x0,y0,x1,y1] ) )
    lines = np.array(lines)
    nx.draw(B,pos=nodepos,
    	with_labels=True
    	#edge_color=np.random.random(nx.number_of_edges(B)),
    	#edge_cmap=plt.get_cmap('Blues'), 
    	#node_color=np.random.random(nx.number_of_nodes(B))
    	#cmap=plt.get_cmap('Reds')
    	)
    plt.show()
    print fi.num_crosses(lines)
    

def drawbipartite2(A):
    from matplotlib.path import Path
    # Create the nx bipartite graph
    SP = np_to_scipy_matrix(A)
    B = nx.bipartite.from_biadjacency_matrix(SP,create_using=None)
    # Compute the nodes positions
    X,Y = nx.bipartite.sets(B)
    pos = dict()
    
    nx.draw(B,pos=pos)
    
    plt.show()

if __name__=='__main__':
    
    filename = sys.argv[1]
    A=np.loadtxt(filename, delimiter=',')
    writedot(A,filename[0:-4])
    writepajek(A,filename[0:-4])
    drawbipartite(A)