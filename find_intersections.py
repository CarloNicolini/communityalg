from __future__ import print_function
import sys
import numpy as np
import math
import random
from random import shuffle
from simanneal import Annealer
#import matplotlib.pyplot as plt

def crosses(lines):
    from shapely.geometry import LineString, Point
    # returns the number of crossings between lines in the list
    points = []
    for i,l1 in enumerate(lines):
        line_shapely1 = LineString([(l1[0],l1[1]),(l1[2],l1[3])])
        for j,l2 in enumerate(lines):
            if j>i:
                line_shapely2 = LineString([(l2[0],l2[1]),(l2[2],l2[3])])
                p = line_shapely1.intersection(line_shapely2)
                
                if not p.is_empty and isinstance(p,Point):
                    points.append(p)
    return points

def num_crosses(lines):
    return len(crosses(lines))


class CrossMinimizer(Annealer):
    def __init__(self, lines,nodes1,nodes2):
        self.lines = lines
        self.nodes1 = nodes1
        self.nodes2 = nodes2
        # state contains the variables of y position of nodes1,nodes2
        # as follows:
        # | y_(n1)_1
        # | y_(n1)_2
        # | y_(n1)_3
        # | y_(n1)_|n1|
        # | .....
        # | y_(n2)_1
        # | y_(n2)_2
        # | y_(n2)_3
        # | y_(n2)_|n2|
        super(CrossMinimizer,self).__init__(nodes1+nodes2)

    def move(self):
        """ Swaps two nodes on the left and on the right """
        r = np.array(xrange(0,len(self.nodes1)-1))
        shuffle(r)
        self.state[range(0,len(self.nodes1)-1)] = self.state[r]
        

    def energy(self):
        L = np.array(self.lines)
        n1 = len(self.nodes1)
        n2 = len(self.nodes2)

        L = np.hstack( (L[self.state[0:n1-1],0:1],L[self.state[n1+1:n2-1],2:3]) )
        return num_crosses(L)


if __name__ == '__main__':
    nrows=10
    A=np.random.randint(-5,5,[nrows,4])
    # Fix the left and right columns
    A[:,0]=0
    A[:,2]=5

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
    init_state
    crossmin = CrossMinimizer(init_state,A,[1,2,3],[4,5])
    crossmin.copy_strategy = "slice"
    state,crosses = crossmin.anneal()    
