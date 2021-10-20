from gurobipy import *
import random
import math
import numpy as np
from queue import PriorityQueue
import time
from SupportingFunctions import *
import collections
import sys
import datetime

"""
Set up for tests
"""

alpha = 50
NHubs = 20
LinkFiltering = False
TripFiltering = False
RunMIP = True
day = "Monday-1.csv"

if len(sys.argv) > 1:
    day = str(sys.argv[2])
    alpha = float(sys.argv[3])
    NHubs = int(sys.argv[4])
    TripFiltering = bool(int(sys.argv[5]))
    LinkFiltering = bool(int(sys.argv[6]))


Stats = collections.defaultdict(int)
Stats['Alpha'] = alpha
Stats['Iterations'] = 0
Stats['LazyCuts'] = 0
Stats['Hubs'] = NHubs
Stats['TripFiltering'] = TripFiltering
Stats['LinkFiltering'] = LinkFiltering
Stats['Day'] = day
alpha = alpha/3600 #convert $/hr to $/sec
"""
Initialize sets and import data
"""

# hubs may be 10 or 20
if NHubs == 10:
    H = ImportHubs('10hubs.dat')
elif NHubs == 20:
    H = ImportHubs('20hubs.dat')
elif NHubs == 30:
    H = ImportHubs('30hubs.dat')
else:
    crash
    
c = 0.00196 # Cost per m
b = 0.0045 # Cost per m
S = 30  # wait time in seconds
n = 32

EPS = 0.0001

# Distances always uses disttime
N,Tau,Gamma,Beta = BuildCosts('disttime.csv',c,b,S,n,alpha,H)

# Day could be one of many files
T = ImportTrips(day,Tau)

print(len(T))

print('data loaded')

# Construct D^r:   
D = {}
for Or,De in T:
    D[Or,De] = [(Or,De)] # direct taxi trip from or to de
    for h in H: # for each hub
        if Or!=h and De != h:
            D[Or,De].append((Or,h))
            D[Or,De].append((h,De))
    
""" 
Trip Filtering & Building bundles of Legs
"""
if TripFiltering:
    SHORTEST, LEGS, removed = TripsAndLegs(T,H,Tau,Gamma,Beta,D)

    print('bundles of legs:',len(LEGS))

    print(len(D),len(removed),len(D)+len(removed))

"""
Link Filtering
"""

if LinkFiltering:
    Links(D, H, Tau, Gamma)

"""
Set up MIP variables and objective
"""

start = time.time()
if RunMIP :
    MIP = Model("MIP")

    X = {(i,j,r): MIP.addVar() for r,dr in D.items() for (i,j) in dr}
    Y = {(h,l,r): MIP.addVar() for r in D for (h,l) in Beta}
    Z = {(h,l): MIP.addVar(vtype=GRB.BINARY) for (h,l) in Beta}

    MIP.setObjective(quicksum(Tau[i,j]*X[i,j,r]*T[r] for (i,j,r) in X)
                    +quicksum(Gamma[h,l]*Y[h,l,r]*T[r] for (h,l,r) in Y)
                    +quicksum(Beta[h,l]*Z[h,l] for (h,l) in Z)
                    ,GRB.MINIMIZE)

    WeakConnectivity = {h:MIP.addConstr(
            quicksum(Z[h,l] for l in H if l != h) ==
            quicksum(Z[l,h] for l in H if l != h)) 
        for h in H}

    OnlyUseOpenBusArcs = {(h,l,r): MIP.addConstr(Y[h,l,r] <= Z[h,l]) for (h,l,r) in Y}

    # Taxi Arc Flow Conservation Constraints (for when i IS NOT in H)
    OrConserve = {(Or,De): MIP.addConstr(
            quicksum(X[i,j,(Or,De)] for (i,j) in dr if i == Or) +
            quicksum(Y[Or,h,(Or,De)]-Y[h,Or,(Or,De)] 
                for h in H if Or in H and Or!=h)==1) 
            for (Or,De),dr in D.items()}

    DeConserve = {(Or,De): MIP.addConstr(
            quicksum(X[i,j,(Or,De)] for (i,j) in dr if j == De) -
            quicksum(Y[De,h,(Or,De)]-Y[h,De,(Or,De)] 
                for h in H if De in H and De!=h)==1)
            for (Or,De),dr in D.items()}


    HubConserve = {((Or,De),i): MIP.addConstr(
            (X[Or,i,(Or,De)] if (Or,i) in dr else 0) 
            -(X[i,De,(Or,De)] if (i,De) in dr else 0)+
            quicksum(Y[h,i,(Or,De)]-Y[i,h,(Or,De)] for h in H if i!=h) ==0) 
            for i in H for (Or,De),dr in D.items() if i!=Or and i!=De}

    MIP.setParam('Threads',4)
    MIP.setParam('TimeLimit',1200.0)
    
loadtime = time.time() - start

Stats['LoadTime'] = loadtime

if RunMIP:
    MIP.optimize()
    Stats['ObjVal'] = MIP.objVal
    
Stats['RunTime'] = time.time() - start

#for z in Z:
#    if Z[z].x > 0.9:
#        print(Z[z].x, z)

if len(sys.argv) > 1:
    with open(str(sys.argv[1]), "a") as f:
        f.write(sys.argv[0]+','+str(datetime.datetime.now())+','+
                str(sorted([(k,v) for (k,v) in zip(Stats.keys(),Stats.values())]))+'\n')  
