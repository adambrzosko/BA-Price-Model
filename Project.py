#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 23:54:34 2021

@author: adam
"""
#import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pickle
import collections
from scipy import stats
#import random
#import os 
outputdir = '/home/adam/Desktop/work/3rd Year/Networks/Project' 

#1. think about proper testing - adjacency matrix for small network, weight etc?
#2. derive the degree distribution by hand for all three phases
#3. what statistical tests to use? KS - too few values? chisquared - not random
#4. for part 4 you want the largest m to show discrepancy?
#5. errors and error bars
#6. how to explain the logbinning motivation and the stats difference

StartingSize = 64
network = [[] for i in range(StartingSize)]
for i in range(StartingSize):
    att = np.random.choice(64, np.random.choice(64), replace = False)
    for k in att:
        if i not in network[k]:
            network[k].append(i)
            network[i].append(k)

network = [[1,2,4],[0,2,3],[0,1,3],[1,2],[0]]           #initial network
weights = []                            #how many edges each vertex has
X = len(network)                        #size of initial network
T = 0                                   #time
#attachments = []
for i in range(X):
    weights.append(len(network[i]))
#    for k in range(len(network[i])):
#        attachments.append(i)
Weight = sum(weights)                   #sum of all edges 


P = []
for i in range(X):
    P.append(weights[i]/Weight)         #probabilities of assignment (preferential)


N = 10000 - X                              #final number of vertices
S = 3                                   #how many edges to add on each run

def AddNode():
    global T
    global Weight
    T += 1
    network.append([])
    weights.append(0)

def AddEdgesPref(m = 3):
    global Weight
    for i in range(T+X-1):                         #recalculate probabilities
        P[i] = weights[i]/Weight   
                              
    #preferentially choose which old ones to connect with the new one
    vertices = np.random.choice(T+X-1, m, replace = False, p=P)  
#    vertices = []
#    while len(vertices) < m:
#        a = np.random.choice(attachments)
#        if a not in vertices:
#            vertices.append(a)
#        else:
#            pass
    for k in range(m):                 #make connections and weight ammendments                                           
        network[vertices[k]].append(X+T-1)
        network[-1].append(k)
        weights[vertices[k]] += 1
        weights[-1] += 1
#        attachments.append(vertices[k])
#        attachments.append(X+T-1)
    Weight += 2*m
    P.append(0)             #append any new P so that the array is right size
                            #(gets recalculated on next iteration anyway)
                            
def DrivePref(s = 3):          #run with s as number of new edges per node added
    global Weight
    AddNode()
    AddEdgesPref(m = s)
    
def AddEdgesRandom(m = 3):
    global Weight  
                              
    #randomly choose which old ones to connect with the new one
    vertices = np.random.choice(T+X-1, m, replace = False)  

    for k in range(m):                 #make connections and weight ammendments                                           
        network[vertices[k]].append(X+T-1)
        network[-1].append(k)
        weights[vertices[k]] += 1
        weights[-1] += 1
    Weight += 2*m
    P.append(0)
                            
def DriveRandom(s = 3):          #run with s as number of new edges per node added
    global Weight
    AddNode()
    AddEdgesRandom(m = s)

def DriveMixed(s = 3, q=2/3):
    '''
    Attach new edge preferentially (prob q) or randomly (prob 1-q)
    '''
    global Weight
    AddNode()
    if np.random.choice([0,1], p=[q,1-q]) == 1:
        AddEdgesRandom(m = s)
    else:
        AddEdgesPref(m = s)

def logbin(data, scale = 1., zeros = False):
    """
    Taken from:
    Max Falkenberg McGillivray
    mff113@ic.ac.uk
    2019 Complexity & Networks course
    """
    if scale < 1:
        raise ValueError('Function requires scale >= 1.')
    count = np.bincount(data)
    tot = np.sum(count)
    smax = np.max(data)
    if scale > 1:
        jmax = np.ceil(np.log(smax)/np.log(scale))
        if zeros:
            binedges = scale ** np.arange(jmax + 1)
            binedges[0] = 0
        else:
            binedges = scale ** np.arange(1,jmax + 1)
            # count = count[1:]
        binedges = np.unique(binedges.astype('uint64'))
        x = (binedges[:-1] * (binedges[1:]-1)) ** 0.5
        y = np.zeros_like(x)
        count = count.astype('float')
        for i in range(len(y)):
            y[i] = np.sum(count[binedges[i]:binedges[i+1]]/(binedges[i+1] - binedges[i]))
    else:
        x = np.nonzero(count)[0]
        y = count[count != 0].astype('float')
        if zeros != True and x[0] == 0:
            x = x[1:]
            y = y[1:]
    y /= tot
    x = x[y!=0]
    y = y[y!=0]
    return x,y

def DegDist(x, m=1):
    p = 2*(m)*(m+1)/(x*(x+1)*(x+2))
    return p

#%%
'''
task 1.3 data
'''
for k in range(1,5):
    network = [[1,2,4],[0,2,3],[0,1,3],[1,2],[0]]           #initial network
    weights = []                            #how many edges each vertex has
    X = len(network)                        #size of initial network
    T = 0                                   #time
    for i in range(X):
        weights.append(len(network[i]))
    Weight = sum(weights)                   #sum of all edges 
    P = []
    for i in range(X):
        P.append(weights[i]/Weight)         #probabilities of assignment (preferential)


    for i in range(1000000):
        DrivePref(k)
    with open(outputdir+'/deg1M-M{}'.format(k), 'wb') as f:
        pickle.dump(weights,f)
    with open(outputdir+'/net1M-M{}'.format(k), 'wb') as f:
        pickle.dump(network,f)      
        
#%%
'''
task 1.3
'''

C = globals()
for k in range(1,5):
    with open(outputdir+'/deg100k-M{}'.format(k), 'rb') as f:
         C['C{}'.format(k)] = sorted(collections.Counter(pickle.load(f)).items())
    x, y = zip(*C['C{}'.format(k)])
    x = np.array(x)
    total = sum(y)
    y = np.array([i/total for i in y])
    plt.plot(x,y, label = 'M = {}'.format(k))
    
    test_stat = stats.ks_2samp(y, DegDist(x, m=k))
    print(test_stat) #too few data points for KS??
    slope, intercept, r_value, p_value, std_err = stats.linregress(y, DegDist(x, m=k))
    print(slope, intercept, r_value, p_value, std_err)
    
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()
   
#%%
'''
logbinned
'''
X = globals()
Y = globals()
for k in range(1,5):
    with open(outputdir+'/deg100k-M{}'.format(k), 'rb') as f:
         X['X{}'.format(k)], Y['Y{}'.format(k)] = logbin(pickle.load(f), 1.3)
         plt.plot(X['X{}'.format(k)], Y['Y{}'.format(k)], '-x', label = 'M = {}'.format(k))
         
         #plt.plot(X['X{}'.format(k)], DegDist(X['X{}'.format(k)], m=k), '-', label = 'D = {}'.format(k))
         test_stat = stats.ks_2samp(Y['Y{}'.format(k)], DegDist(X['X{}'.format(k)], m=k))
         print(test_stat) #too few data points for KS??
         slope, intercept, r_value, p_value, std_err = stats.linregress(Y['Y{}'.format(k)], DegDist(X['X{}'.format(k)], m=k))
         print(slope, intercept, r_value, p_value, std_err)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()

#%%
'''
1.4 data
we go to highest m possible, m=5
'''
for j in range(10,18):
    network = [[1,2,4],[0,2,3],[0,1,3],[1,2],[0]]           #initial network
    weights = []                            #how many edges each vertex has
    X = len(network)                        #size of initial network
    T = 0                                   #time
    for i in range(X):
        weights.append(len(network[i]))
    Weight = sum(weights)                   #sum of all edges 
    P = []
    for i in range(X):
        P.append(weights[i]/Weight)         #probabilities of assignment (preferential)
    
    for k in range(2**j): #do like power of 2 or something
        DrivePref(s=5)
    with open(outputdir+'/degM5-N2e{}'.format(j), 'wb') as f:
        pickle.dump(weights,f)
    with open(outputdir+'/netM5-N2e{}'.format(j), 'wb') as f:
        pickle.dump(network,f)     




