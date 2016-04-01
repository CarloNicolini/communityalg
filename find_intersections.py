#from __future__ import print_function
import sys
import numpy as np
import math
import random
from random import shuffle
from simanneal import Annealer
import networkx as nx
from scipy.sparse import csr_matrix as np_to_scipy_matrix

eps = np.finfo(float).eps
#import matplotlib.pyplot as plt


class CrossMinimizer(Annealer):
    def __init__(self, state, graph):
        self.graph = graph # a networkx Graph object to represent a bipartite graph
        self.state = state
        self.nodesleft, self.nodesright = nx.bipartite.sets(self.graph)
        self.nA = len(self.nodesleft)
        self.nB = len(self.nodesright)
        
        self.nodepos={}
        self.nodepos.update( [n,(-1,max(self.nA,self.nB)-i) ] for i,n in enumerate(self.nodesleft))
        self.nodepos.update( [n,(1,max(self.nA,self.nB)-i) ] for i,n in enumerate(self.nodesright))    
        # Generate the set of lines for all the edges
        self.lines = []
        for e in graph.edges():
            n0 = self.state[e[0]]
            n1 = self.state[e[1]]
            x0 = self.nodepos[n0][0]
            y0 = self.nodepos[n0][1]
            x1 = self.nodepos[n1][0]
            y1 = self.nodepos[n1][1]
            self.lines.append( np.array([x0,y0,x1,y1] ) )
        self.lines = np.array(self.lines)
        super(CrossMinimizer,self).__init__(state)

    def crosses(self,lines):
        from shapely.geometry import LineString, Point
        # returns the number of crossings between lines in the list
        points = []
        for i,l1 in enumerate(lines):
            line_shapely1 = LineString([(l1[0],l1[1]),(l1[2],l1[3])])
            for j,l2 in enumerate(lines):
                # avoid points that share the same extremal point
                if j>i and (l1[1]==l2[1] ) : # qui mettere le linee non dello stesso vertice
                    line_shapely2 = LineString([(l2[0],l2[1]),(l2[2],l2[3])])
                    p = line_shapely1.intersection(line_shapely2)
                    if not p.is_empty and isinstance(p,Point):
                        points.append(p)
        return points

    def crosses2(self):
        from shapely.geometry import LineString, Point
        

    def num_crosses(self,lines):
        return len(self.crosses(lines))

    def move(self):
        """ Swaps two nodes on the left and on the right """
        # Choose two random nodes in the nodesetA to swap
        nodes_swap_left = np.random.random_integers(0,self.nA-1,2)
        nodes_swap_left=[0,3]
        print "Swapping %d (%s) with %d (%s)" % (nodes_swap_left[0],self.nodepos[nodes_swap_left[0]],nodes_swap_left[1],self.nodepos[nodes_swap_left[1]])
        # Do the swap
        self.state[nodes_swap_left[0]],self.state[nodes_swap_left[1]] = self.state[nodes_swap_left[1]],self.state[nodes_swap_left[0]]
        
        # Choose two random nodes in the nodesetA to swap
        #nodes_swap_right = np.random.random_integers(self.nA,self.nA+self.nB-1,2)
        # Do the swap
        # Update nodepos
        self.nodepos[nodes_swap_left[0]],self.nodepos[nodes_swap_left[1]] =         self.nodepos[nodes_swap_left[1]],self.nodepos[nodes_swap_left[0]]
        #self.state[nodes_swap_right[0]],self.state[nodes_swap_right[1]] = self.state[nodes_swap_right[1]],self.state[nodes_swap_right[0]]
        # Update the self.line structure
        self.lines = []
        for e in self.graph.edges():
            n0 = self.state[e[0]]
            n1 = self.state[e[1]]
            x0 = self.nodepos[n0][0]
            y0 = self.nodepos[n0][1]
            x1 = self.nodepos[n1][0]
            y1 = self.nodepos[n1][1]
            self.lines.append( np.array([x0,y0,x1,y1] ) )
        self.lines = np.array(self.lines)

    def energy(self):
        return self.num_crosses(self.lines)

    def plot(self):
        import matplotlib.pyplot as plt
        print self.nodepos
        nx.draw(self.graph,self.nodepos,with_labels=True)
        plt.show()

    def info(self):
        print("Current crossings = %d") % self.energy()
        for n,xy in self.nodepos.iteritems():
            print("Node %d x=%d y=%d") % (n,xy[0],xy[1])
        for r in self.lines:
            print r


if __name__ == '__main__':
    filename = sys.argv[1]
    A = np.loadtxt(filename,delimiter=',')
    G=nx.bipartite.from_biadjacency_matrix(np_to_scipy_matrix(A))

    #G = nx.bipartite.from_biadjacency_matrix(np.loadtxt(filename, delimiter=','))
    # plt.subplots(2,sharex=True)
    
    # for i in range(0,nrows):
    #     plt.plot([A[i,0],A[i,2]],[A[i,1],A[i,3]])

    # points = crosses(A)
    # plt.plot([p.x for p in points],[p.y for p in points],'ko')

    # plt.draw()
    # plt.show()

    # init_state = range(0,nrows)

    # points = crosses(A)

    # #plt.draw()
    # #plt.show()
    # #sys.exit(0)
    
    crossmin = CrossMinimizer(range(0,7),G)
    crossmin.copy_strategy = "slice"
    #X = crossmin.auto(minutes=1,steps=100)
    #crossmin.plot()
    crossmin.info()
    crossmin.move()

    print "========"
    crossmin.info()
    crossmin.plot()
    # print crossmin.state
    # print crossmin.lines
    # print "#crossings = ",crossmin.energy()
    # crossmin.plot()
    #state, ncrossings = crossmin.anneal()