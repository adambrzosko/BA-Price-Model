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


StartingSize = 64
#network = [[] for i in range(StartingSize)]
#for i in range(StartingSize):
#    att = np.random.choice(StartingSize, np.random.choice(np.arange(1,StartingSize)), replace = False)
#    for k in att:
#        if i not in network[k] and i!=k:
#            network[k].append(i)
#            network[i].append(k)
            
network = [[] for i in range(StartingSize)]
for i in range(StartingSize):
    for k in range(StartingSize):
        if i!=k:    
            network[i].append(k)

#network = [[1,2,4],[0,2,3],[0,1,3],[1,2],[0]]           #initial network
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


def logbin(data, scale = 1., MaxDeg = 0, zeros = True):
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
    if MaxDeg > 0:
        smax = MaxDeg
    else:
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
for it in range(1,11):
    for r in [2, 4, 8, 16, 32]:
        network = [[] for i in range(r)]
        for i in range(r):
            for k in range(r):
                if i!=k:    
                    network[i].append(k)

        weights = []                            #how many edges each vertex has
        X = len(network)                        #size of initial network
        T = 0                                   #time
        for i in range(X):
            weights.append(len(network[i]))
        Weight = sum(weights)                   #sum of all edges 
        P = []
        for i in range(X):
            P.append(weights[i]/Weight)         #probabilities of assignment (preferential)


        for i in range(100000):
            DrivePref(r)
        with open(outputdir+'/1.3Final/deg100k-M{}Run{}'.format(r,it), 'wb') as f:
            pickle.dump(weights,f)
        with open(outputdir+'/1.3Final/net100k-M{}Run{}'.format(r,it), 'wb') as f:
            pickle.dump(network,f)      
        
#%%
'''
task 1.3
'''

C = globals()
for k in range(1,6):
    with open(outputdir+'/deg100k-M{}Run4'.format(k), 'rb') as f:
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
Data = {}
for k in range(1,6):
    L = 2**k
    Maxes = []
    for it in range(2,10):
        with open(outputdir+'/1.3RandomInitial100k/deg100k-M{}Run{}'.format(L,it), 'rb') as f:
            Data['M{}Run{}'.format(L,it)] = pickle.load(f)
            Maxes.append(np.max(Data['M{}Run{}'.format(L,it)]))
    Data['KmaxM{}'.format(L)] = np.max(Maxes)
    
for k in range(1,6):
    L = 2**k
    Data['X{}'.format(L)] = []
    Data['Y{}'.format(L)] = []
    Biggest = 0
    for it in range(2,10):
        Data['X{}Run{}'.format(L,it)], Data['Y{}Run{}'.format(L,it)] = logbin(Data['M{}Run{}'.format(L,it)], 1.1, Data['KmaxM{}'.format(L)])
       
        for x in Data['X{}Run{}'.format(L,it)]:
            if x not in Data['X{}'.format(L)]:
                Data['X{}'.format(L)].append(x)
                Data['Y{}'.format(L)].append(Data['Y{}Run{}'.format(L,it)][np.ndarray.tolist(Data['X{}Run{}'.format(L,it)]).index(x)])
            else:
                Data['Y{}'.format(L)][Data['X{}'.format(L)].index(x)] += Data['Y{}Run{}'.format(L,it)][np.ndarray.tolist(Data['X{}Run{}'.format(L,it)]).index(x)]
                

for k in range(1,6):
    L = 2**k
    Data['Y{}'.format(L)] = [x/9 for x in Data['Y{}'.format(L)]]
    plt.plot(Data['X{}'.format(L)],Data['Y{}'.format(L)], '+', label = 'M = {}'.format(L))

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
    StartingSize = 32
    network = [[] for i in range(StartingSize)]
    for i in range(StartingSize):
        att = np.random.choice(StartingSize, np.random.choice(np.arange(1,StartingSize)), replace = False)
        for k in att:
            if i not in network[k] and i!=k:
                network[k].append(i)
                network[i].append(k)
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
        DriveRandom(s=2)
    with open(outputdir+'/MixdegM2-N2e{}'.format(j), 'wb') as f:
        pickle.dump(weights,f)
    with open(outputdir+'/MixdnetM2-N2e{}'.format(j), 'wb') as f:
        pickle.dump(network,f)     

