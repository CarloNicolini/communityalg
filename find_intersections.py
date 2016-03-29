from __future__ import print_function
import sys
import numpy as np
import math
import random
from random import shuffle
from simanneal import Annealer
import matplotlib.pyplot as plt

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

def num_crosses(lines,ordering):
    return len(crosses(lines[ordering,:]))


class CrossMinimizer(Annealer):
    def __init__(self, state, lines):
        self.lines = lines
        super(CrossMinimizer,self).__init__(state)

    def move(self):
        """ Swaps two nodes on the left and on the right """
        a = random.randint(0,len(self.state)-1)
        b = random.randint(0,len(self.state)-1)
        #print((a,b))
        # Do the actual swap
        self.state[a],self.state[b] = self.state[b],self.state[a]

    def energy(self):
        return num_crosses(self.lines,self.state)


if __name__ == '__main__':
    nrows=5
    A=np.random.randint(-5,5,[nrows,4])
    # Fix the left and right columns
    A[:,0]=0
    A[:,2]=5

    plt.subplots(2,sharex=True)
    
    for i in range(0,nrows):
        plt.plot([A[i,0],A[i,2]],[A[i,1],A[i,3]])

    points = crosses(A)
    plt.plot([p.x for p in points],[p.y for p in points],'ko')

    plt.draw()
    plt.show()

    init_state = range(0,nrows)
    shuffle(init_state)
    print(init_state)
    A = A[init_state,:]

    for i in range(0,nrows):
        plt.plot([A[i,0],A[i,2]],[A[i,1],A[i,3]])

    points = crosses(A)

    for i in range(0,len(points)):
        plt.plot(points[i].x,points[i].y,'ko')

    plt.draw()
    plt.show()
    sys.exit(0)

    crossmin = CrossMinimizer(init_state,A)
    crossmin.copy_strategy = "slice"
    state,crosses = crossmin.anneal()    
