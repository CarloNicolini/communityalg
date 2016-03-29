#!/usr/bin/python

import sys
import numpy as np
import networkx as nx
from scipy.sparse import csr_matrix as np_to_scipy_matrix
import matplotlib as mpl
import matplotlib.pyplot as plt

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

    # Compute the nodes positions
    X,Y = nx.bipartite.sets(B)
    n1,n2 = len(X),len(Y)
    nodepos = dict()
    nodepos.update( [n,(-1,i) ] for i,n in enumerate(X))
    nodepos.update( [n,(1,i) ] for i,n in enumerate(Y))
    print nodepos
    """
    nx.draw(B,pos=nodepos,
    	with_labels=True
    	#edge_color=np.random.random(nx.number_of_edges(B)),
    	#edge_cmap=plt.get_cmap('Blues'), 
    	#node_color=np.random.random(nx.number_of_nodes(B))
    	#cmap=plt.get_cmap('Reds')
    	)
    """
    #import matplotlib as mpl
    #import matplotlib.pyplot as plt
    #plt.show()
    edges = []
    for node1,p in nodepos.iteritems():
        for node2,p in nodepos.iteritems():
            if B.is_edge(node1,node2):
                edges.append([])
                
    all_edges_pos = np.zeros(m,4)
    print X
    all_edges_pos[X]

def perp( a ) :
    b = empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

# line segment a given by endpoints a1, a2
# line segment b given by endpoints b1, b2
# return 
def seg_intersect(a1,a2, b1,b2):
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = perp(da)
    denom = dot( dap, db)
    num = dot( dap, dp )
    return (num / denom.astype(float))*db + b1
    

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