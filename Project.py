#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 23:54:34 2021

@author: adam
"""
#import networkx as nx
import numpy as np
#import matplotlib.pyplot as plt
#import random
#import os 
outputdir = '/home/adam/Desktop/work/3rd Year/Networks' 

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


N = 1000 - X                              #final number of vertices
S = 3                                   #how many edges to add on each run

def AddNode():
    global T
    global Weight
    T += 1
    network.append([])
    weights.append(0)

def AddEdges(m = 3):
    global Weight
    for i in range(T+X-1):                         #recalculate probabilities
        P[i] = weights[i]/Weight   
                              
    #preferentially choose which old ones to connect with the new one
    vertices = np.random.choice(T+X-1, m, replace = False, p=P)  
    
    for k in range(m):                 #make connections and weight ammendments                                           
        network[vertices[k]].append(X+T-1)
        network[-1].append(k)
        weights[vertices[k]] += 1
        weights[-1] += 1
    Weight += 2*m
    P.append(0)             #append any new P so that the array is right size
                            #(gets recalculated on next iteration anyway)
                            
def Drive(s = 3):          #run with s as number of new edges per node added
    global Weight
    AddNode()
    AddEdges(m = s)

#%%
for i in range(N):  
    Drive(S)
    
