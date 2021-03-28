#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 23:54:34 2021

@author: adam
"""
#import networkx as nx
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt
import pickle
from scipy import stats
from scipy.optimize import curve_fit as cf
#import random
#import os 
outputdir = '/home/adam/Desktop/work/3rd Year/Networks/ProjectData' 
graphdir = '/home/adam/Desktop/work/3rd Year/Networks/Report'


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
attachments = []
for i in range(X):
    weights.append(len(network[i]))
    for k in range(len(network[i])):
        attachments.append(i)
Weight = sum(weights)                   #sum of all edges 

def AddNode():
    global T
    global Weight
    T += 1
    network.append([])
    weights.append(0)
    
def AddEdgesPref(m = 3):
    global Weight 
                              
    vertices = []
    while len(vertices) < m:
        a = attachments[random.randint(len(attachments))]
        if a not in vertices:
            vertices.append(a)
        else:
            pass
        
    for k in range(m):                 #make connections and weight ammendments                                           
        network[vertices[k]].append(X+T-1)
        network[-1].append(vertices[k])
        weights[vertices[k]] += 1
        weights[-1] += 1
        attachments.append(vertices[k])
        attachments.append(X+T-1)
    Weight += 2*m

                            
def DrivePref(s = 3):          #run with s as number of new edges per node added
    global Weight
    AddNode()
    AddEdgesPref(m = s)
    
def AddEdgesRandom(m = 3):
    global Weight  
                              
    #randomly choose which old ones to connect with the new one
    vertices = []
    while len(vertices) < m:
        a = random.randint(T+X-1)
        if a not in vertices:
            vertices.append(a)
        else:
            pass

    for k in range(m):                 #make connections and weight ammendments                                           
        network[vertices[k]].append(X+T-1)
        network[-1].append(vertices[k])
        weights[vertices[k]] += 1
        weights[-1] += 1
        attachments.append(vertices[k])
        attachments.append(X+T-1)
    Weight += 2*m
                            
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
    if random.choice([0,1], p=[q,1-q]) == 1:
        AddEdgesRandom(m = s)
    else:
        AddEdgesPref(m = s)


def logbin(data, scale = 1., MaxDeg = 0, MinDeg = 2, zeros = True):
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
        jmin = np.floor(np.log(MinDeg)/np.log(scale))
        if zeros:
            binedges = scale ** np.arange(jmin, jmax + 1)  #was just jmax+1
            binedges[0] = MinDeg                #was 0
        else:
            binedges = scale ** np.arange(1,jmax + 1)
            # count = count[1:]
        binedges = np.unique(binedges.astype('uint64'))
        x = (binedges[:-1] * (binedges[1:]-1)) ** 0.5
        y = np.zeros_like(x)
        stdev = np.zeros_like(x)
        number = np.zeros_like(x)
        count = count.astype('float')
        for i in range(len(y)):
            y[i] = np.sum(count[binedges[i]:binedges[i+1]]/(binedges[i+1] - binedges[i]))
            number[i] = np.sum(count[binedges[i]:binedges[i+1]]) #number of avalnches in bin i
#        y /= tot
#        for i in range(len(y)):
#            probs = []
#            for j in range(binedges[i]-1,binedges[i+1]-1):
#                probs.append((y[i]-(count[j]/tot))**2)
#            stdev[i] = np.sqrt(sum(probs)/(binedges[i+1] - binedges[i]))
    else:
        x = np.nonzero(count)[0]
        y = count[count != 0].astype('float')
        if zeros != True and x[0] == 0:
            x = x[1:]
            y = y[1:]
#        y /= tot
#        for i in range(len(y)):
#            for j in range(binedges[i],binedges[i+1]):
#                stdev[i] = np.sum((y[i]-(count[j]/tot))**2)/(binedges[i+1] - binedges[i])
    y /= tot
    x = x[y!=0]
    y = y[y!=0]
#    stdev = stdev[stdev!=0]
    return x,y,number

def DegDist(x, m=1):
    p = 2*(m)*(m+1)/(x*(x+1)*(x+2))
    return p

def DegDistRand(x,m=1):
    #p = (m**(x-m))/((1+m)**(1+x-m))
    p = ((m/(1+m))**(x-m))*(1/(1+m))
    return p

def DegDistMix12(x, m=1):
    A = 12*m*(3*m+1)*(3*m+2)*(3*m+3)
    B = (k+2*m)*(k+2*m+1)(k+2*m+2)(k+2*m+3)(k+2*m+4)
    p = A/B
    return p

def DegDistMix23(x, m=1):
    A = 6*m*(2*m+1)*(2*m+2)
    B = (k+m)*(k+m+1)(k+m+2)(k+m+3)
    p = A/B
    return p

#%%
'''
task 1.3 data
'''
for it in range(101,201):
    for r in [2, 4, 8, 16, 32, 64]:
        network = [[] for i in range(r)]
        for i in range(r):
            for k in range(r):
                if i!=k:    
                    network[i].append(k)
        weights = []                            
        X = len(network)                        
        T = 0                                   
        attachments = []
        for i in range(X):
            weights.append(len(network[i]))
            for k in range(len(network[i])):
                attachments.append(i)
        Weight = sum(weights)                   

        for i in range(100000):
            DriveMixed(r,0.5)
        with open(outputdir+'/3.3Final12/deg100k-M{}Run{}'.format(r,it), 'wb') as f:
            pickle.dump(weights,f)
#        with open(outputdir+'/1.3Final/net100k-M{}Run{}'.format(r,it), 'wb') as f:
#            pickle.dump(network,f)      
        
   
#%%
'''
1.3 analysis
'''
colors = ['magenta', 'firebrick', 'royalblue', 'orangered', 'black', 'forestgreen', 'navy']
Data = {}    #dictionary for all the data
for k in range(1,7):
    L = 2**k
    Maxes = []    #list of max degrees in every run
#    Data['All{}'.format(L)] = [] #a list for all degrees from all runs
    for it in range(1,201):
        with open(outputdir+'/1.3Final/deg100k-M{}Run{}'.format(L,it), 'rb') as f:
            Data['M{}Run{}'.format(L,it)] = pickle.load(f)   #load all runs with different ms
            Maxes.append(np.max(Data['M{}Run{}'.format(L,it)]))  #find max degree for every run and every m
#            Data['All{}'.format(L)] += Data['M{}Run{}'.format(L,it)]
    Data['KmaxM{}'.format(L)] = np.max(Maxes)  #get the ultimate max degree for every m
#    Data['Xall{}'.format(L)], Data['Yall{}'.format(L)], Data['Numall{}'.format(L)] = logbin(Data['All{}'.format(L)], 1.3, Data['KmaxM{}'.format(L)], L)
for k in range(1,7):
    L = 2**k
    Data['X{}'.format(L)] = []
    Data['Y{}'.format(L)] = []
    Data['Number{}'.format(L)] = []
#    Data['Dev{}'.format(L)] = []
#    Data['YallStdev{}'.format(L)] = [0 for x in Data['Xall{}'.format(L)]]
#    Data['YallReads{}'.format(L)] = [0 for x in Data['Xall{}'.format(L)]]
    for it in range(1,201): #logbin the data for each run and each m
        Data['X{}Run{}'.format(L,it)], Data['Y{}Run{}'.format(L,it)], Data['Number{}Run{}'.format(L,it)] = logbin(Data['M{}Run{}'.format(L,it)], 1.3, Data['KmaxM{}'.format(L)], L)
       
        for x in Data['X{}Run{}'.format(L,it)]:  #run through all runs of different m 
            if x not in Data['X{}'.format(L)]:   #make sure the master array for each m has all the bins
                Data['X{}'.format(L)].append(x)  #append their probabilities below if the master doesn't have that entry
                Data['Y{}'.format(L)].append(Data['Y{}Run{}'.format(L,it)][np.ndarray.tolist(Data['X{}Run{}'.format(L,it)]).index(x)])
                Data['Number{}'.format(L)].append(Data['Number{}Run{}'.format(L,it)][np.ndarray.tolist(Data['X{}Run{}'.format(L,it)]).index(x)])
            else: #increase the probabilities in the master by the probs of the run if the bins already are there
                Data['Y{}'.format(L)][Data['X{}'.format(L)].index(x)] += Data['Y{}Run{}'.format(L,it)][np.ndarray.tolist(Data['X{}Run{}'.format(L,it)]).index(x)]
                Data['Number{}'.format(L)][Data['X{}'.format(L)].index(x)] += Data['Number{}Run{}'.format(L,it)][np.ndarray.tolist(Data['X{}Run{}'.format(L,it)]).index(x)]
#                Data['Dev{}'.format(L)][Data['X{}'.format(L)].index(x)] += Data['Dev{}Run{}'.format(L,it)][np.ndarray.tolist(Data['X{}Run{}'.format(L,it)]).index(x)]
#        for stdev in Data['X{}Run{}'.format(L,it)]:
#            Data['YallStdev{}'.format(L)][np.ndarray.tolist(Data['Xall{}'.format(L)]).index(stdev)] += (Data['Y{}Run{}'.format(L,it)][np.ndarray.tolist(Data['X{}Run{}'.format(L,it)]).index(stdev)] - Data['Yall{}'.format(L)][np.ndarray.tolist(Data['Xall{}'.format(L)]).index(stdev)])**2
#            Data['YallReads{}'.format(L)][np.ndarray.tolist(Data['Xall{}'.format(L)]).index(stdev)] += 1
#    Data['Numpy{}'.format(L)] = [[] for i in Data['X{}'.format(L)]]
    Data['Y{}'.format(L)] = [x/200 for x in Data['Y{}'.format(L)]] #divided by the number of runs
    Data['Y{}Var'.format(L)] = [0 for x in Data['X{}'.format(L)]] #for error calculation
    Data['Y{}Read'.format(L)] = [0 for x in Data['X{}'.format(L)]] #also for errors
    for it in range(1,201):
        for x in Data['X{}Run{}'.format(L,it)]: #(below) append the variances in the right places from all the runs
            Data['Y{}Var'.format(L)][Data['X{}'.format(L)].index(x)] += (Data['Y{}Run{}'.format(L,it)][np.ndarray.tolist(Data['X{}Run{}'.format(L,it)]).index(x)] - Data['Y{}'.format(L)][Data['X{}'.format(L)].index(x)])**2
            Data['Y{}Read'.format(L)][Data['X{}'.format(L)].index(x)] += 1  #if variance appended, increase the count N
#            Data['Numpy{}'.format(L)][Data['X{}'.format(L)].index(x)].append(Data['Y{}Run{}'.format(L,it)][np.ndarray.tolist(Data['X{}Run{}'.format(L,it)]).index(x)])
    Data['Y{}Stdev'.format(L)] = [np.sqrt(Data['Y{}Var'.format(L)][i]/Data['Y{}Read'.format(L)][i]) for i in range(len(Data['Y{}Var'.format(L)]))]
    Data['X{}'.format(L)], Data['Y{}'.format(L)], Data['Y{}Stdev'.format(L)] = zip(*sorted(zip(Data['X{}'.format(L)], Data['Y{}'.format(L)], Data['Y{}Stdev'.format(L)])))
#    Data['Error{}'.format(L)] = [np.std(Data['Numpy{}'.format(L)][i]) for i in range(len(Data['Y{}Var'.format(L)]))]
#    Data['Y{}Error'.format(L)] = [np.sqrt(Data['YallStdev{}'.format(L)][i]/Data['YallReads{}'.format(L)][i]) for i in range(len(Data['YallStdev{}'.format(L)]))]
#    Data['Y{}Error'.format(L)] = np.array(Data['Y{}Error'.format(L)])
    Data['X{}'.format(L)] = np.array((Data['X{}'.format(L)]))
    Data['Y{}'.format(L)] = np.array((Data['Y{}'.format(L)]))
    Data['Y{}Stdev'.format(L)] = np.array(Data['Y{}Stdev'.format(L)])
for k in range(1,7):
    #base = np.linspace(k**2,60*k,1000) #random
    base = np.linspace(2**k-0.5*k,4800,10000) #preferential
    #base = np.float128(base)
    L = 2**k
    plt.plot(base, DegDist(base,L), '--', color = colors[k-1])
    plt.errorbar(Data['X{}'.format(L)],Data['Y{}'.format(L)], yerr = Data['Y{}Stdev'.format(L)], label = '$m = {}$'.format(L), capsize = 3, color = colors[k-1])
    
plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize = 12, loc = 'lower left')
plt.xlabel(r'$k$', fontsize = 15)
plt.ylabel(r'$p(k)$', fontsize = 15)
plt.savefig(graphdir+'/13Pref', dpi=500)
#%%
'''
1.3 tests
'''
def chisqr_reduced(obs, exp, error, reads):
    chisqr_reduced = 0
    m = max(reads)
    #print(reads)
    for i in range(len(obs)):
        chisqr = ((obs[i]-exp[i])**2)/(error[i])
        chisqr_reduced += chisqr/(m+1-reads[i])
    return chisqr_reduced
for k in range(1,7):
    #base = np.linspace(k**2,60*k,1000) #random
    base = np.linspace(2**k-0.5*k,4800,10000) #preferential
    #base = np.float128(base)
    b = 8-k  #[:-b]
    L = 2**k
    observed = [Data['Y{}'.format(L)][i]*Data['Number{}'.format(L)][i] for i in range(len(Data['Y{}'.format(L)]))]
    expected = [DegDist(Data['X{}'.format(L)][i], m=L)*Data['Number{}'.format(L)][i] for i in range(len(Data['Y{}'.format(L)]))]
    error = [DegDist(Data['Y{}Stdev'.format(L)][i], m=L)*Data['Number{}'.format(L)][i] for i in range(len(Data['Y{}'.format(L)]))]
    test_stat = stats.ks_2samp(Data['Y{}'.format(L)], DegDist(Data['X{}'.format(L)], m=L))
    print(test_stat)
    slope, intercept, r_value, p_value, std_err = stats.linregress(Data['Y{}'.format(L)], DegDist(Data['X{}'.format(L)], m=L))
    print(slope, intercept, r_value, p_value, std_err)
    chi = chisqr_reduced(observed[:-b], expected[:-b], error[:-b], Data['Y{}Read'.format(L)][:-b])
#    chi = stats.chisquare(Data['Y{}'.format(L)], DegDist(Data['X{}'.format(L)], m=L))
#    print('M =', L, 'chi=', chi, 'p=', stats.distributions.chi2.sf(chi,1))
#%%
'''
1.4 data
we go to smalles m possible, m=2
'''
for it in range(101,201):
    for j in range(10,18):
        network = [[] for i in range(2)]
        for i in range(2):
            for k in range(2):
                if i!=k:    
                    network[i].append(k)
        weights = []                            #how many edges each vertex has
        X = len(network)                        #size of initial network
        T = 0                                   #time
        attachments = []
        for i in range(X):
            weights.append(len(network[i]))
            for k in range(len(network[i])):
                attachments.append(i)
        Weight = sum(weights)                   #sum of all edges 
        
    
        for k in range(2**j): #do like power of 2 or something
            DriveMixed(2)
        with open(outputdir+'/3.4Final/degM2-N2e{}Run{}'.format(j,it), 'wb') as f:
            pickle.dump(weights,f)
#        with open(outputdir+'/1.4Final/netM2-N2e{}Run{}'.format(j,it), 'wb') as f:
#            pickle.dump(network,f)     

#%%
'''
1.4 analysis
'''
DataK = {}
for j in range(10,18):
    N = j
    Maxes = []
    for it in range(1,201):
        with open(outputdir+'/1.4Final/degM2-N2e{}Run{}'.format(N,it), 'rb') as f:
            DataK['N{}Run{}'.format(N,it)] = pickle.load(f)
            Maxes.append(np.max(DataK['N{}Run{}'.format(N,it)]))
    DataK['KmaxN{}'.format(N)] = np.max(Maxes)
    DataK['MaxesN{}'.format(N)] = Maxes

for k in range(10,18):
    N = k
    DataK['X{}'.format(N)] = []
    DataK['Y{}'.format(N)] = []
    for it in range(1,201):
        DataK['X{}Run{}'.format(N,it)], DataK['Y{}Run{}'.format(N,it)], DataK['Number{}Run{}'.format(N,it)] = logbin(DataK['N{}Run{}'.format(N,it)], 1.3, DataK['KmaxN{}'.format(N)], 2)
       
        for x in DataK['X{}Run{}'.format(N,it)]:
            if x not in DataK['X{}'.format(N)]:
                DataK['X{}'.format(N)].append(x)
                DataK['Y{}'.format(N)].append(DataK['Y{}Run{}'.format(N,it)][np.ndarray.tolist(DataK['X{}Run{}'.format(N,it)]).index(x)])
            else:
                DataK['Y{}'.format(N)][DataK['X{}'.format(N)].index(x)] += DataK['Y{}Run{}'.format(N,it)][np.ndarray.tolist(DataK['X{}Run{}'.format(N,it)]).index(x)]
    DataK['Y{}'.format(N)] = [x/200 for x in DataK['Y{}'.format(N)]]
    DataK['Y{}Var'.format(N)] = [0 for x in DataK['X{}'.format(N)]]
    DataK['Y{}Read'.format(N)] = [0 for x in DataK['X{}'.format(N)]]
    for it in range(1,201):
        for x in DataK['X{}Run{}'.format(N,it)]:
            DataK['Y{}Var'.format(N)][DataK['X{}'.format(N)].index(x)] +=  (DataK['Y{}Run{}'.format(N,it)][np.ndarray.tolist(DataK['X{}Run{}'.format(N,it)]).index(x)] - DataK['Y{}'.format(N)][DataK['X{}'.format(N)].index(x)])**2
            DataK['Y{}Read'.format(N)][DataK['X{}'.format(N)].index(x)] += 1
    DataK['Y{}Stdev'.format(N)] = [np.sqrt(DataK['Y{}Var'.format(N)][i]/DataK['Y{}Read'.format(N)][i]) for i in range(len(DataK['Y{}Var'.format(N)]))]
    DataK['X{}'.format(N)], DataK['Y{}'.format(N)], Data['Y{}Stdev'.format(N)] = zip(*sorted(zip(DataK['X{}'.format(N)], DataK['Y{}'.format(N)], DataK['Y{}Stdev'.format(N)])))
    DataK['X{}'.format(N)] = np.array((DataK['X{}'.format(N)]))
    DataK['Y{}'.format(N)] = np.array((DataK['Y{}'.format(N)]))
    DataK['Y{}Stdev'.format(N)] = np.array((DataK['Y{}Stdev'.format(N)]))
for k in range(10,18):
    base = np.linspace(k**2,60*k,1000)
    #base = np.linspace(k**2,3000,1000) #preferential
    N = k
    plt.errorbar(DataK['X{}'.format(N)],DataK['Y{}'.format(N)], yerr = DataK['Y{}Stdev'.format(N)], label = 'N = 2^{}'.format(N), capsize = 3)
    #plt.plot(base, DegDistRand(base,L), '-', label = 'Theoretical M = {}'.format(L))
plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize = 12, loc = 'lower left')
plt.xlabel(r'$k$', fontsize = 15)
plt.ylabel(r'$p(k)$', fontsize = 15)
plt.savefig(graphdir+'/14Pref', dpi=500)

#%%
'''
max k
'''
def Fit(N, a, b):
    k = a*(N**b)
    return k

y = []
x = []
Error = []
for k in range(10,18):
    N = k
    AvgMax = np.average(DataK['MaxesN{}'.format(N)])
    Error.append(np.std(DataK['MaxesN{}'.format(N)]))
    y.append(AvgMax)
    x.append(2**N)
    
popt, pcov = cf(Fit, x, y)
print('The slope is', popt[1], 'pm', np.sqrt(pcov[1][1]))
plt.errorbar(x, y, yerr = Error, label = '$k_{1}$', capsize = 3)
plt.plot(x, Fit(x, *popt), '--k', label = 'Two parameter fit')
plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize = 14, loc = 'upper left')
plt.xlabel(r'$k$', fontsize = 15)
plt.ylabel(r'$N$', fontsize = 15)
plt.savefig(graphdir+'/MaxKPref', dpi=500)

#%%
'''
data collapse
'''
for k in range(10,18):
    N = k
    DataK['X{}Scaled'.format(N)] = [DataK['X{}'.format(N)][i]/DataK['KmaxN{}'.format(N)] for i in range(len(DataK['X{}'.format(N)]))]
    DataK['Y{}Scaled'.format(N)] = [DataK['Y{}'.format(N)][i]/DegDist(DataK['X{}'.format(N)][i], m=2) for i in range(len(DataK['X{}'.format(N)]))]
    DataK['Y{}StdevScaled'.format(N)] = [DataK['Y{}Stdev'.format(N)][i]/DegDist(DataK['X{}'.format(N)][i], m=2) for i in range(len(DataK['X{}'.format(N)]))]
    DataK[] #for scaled error
    plt.errorbar(DataK['X{}Scaled'.format(N)],DataK['Y{}Scaled'.format(N)], yerr = , '.', label = 'N = 2^{}'.format(N))
    
plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize = 12, loc = 'lower left')
plt.xlabel(r'$k/k_{1}$', fontsize = 15)
plt.ylabel(r'$p(k)/p_{\infty}(k)$', fontsize = 15)
plt.savefig(graphdir+'/CollapsePref', dpi=500)
