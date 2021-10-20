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
DisaggregatedAnalytic = True
DisaggregatedPSP = False
global AnalyticShortestPath
AnalyticShortestPath = False
ShortestPathPrecuts = True
NHubs = 20
LinkFiltering = False
LinkFiltering2 = False
AnalyticCuts2 = False
Heuristic = False
day = 'Monday-1.csv'

if len(sys.argv) > 1:
    day = str(sys.argv[2])   
    alpha = float(sys.argv[3])
    NHubs = int(sys.argv[4])
    
    DisaggregatedAnalytic = bool(int(sys.argv[5]))
    DisaggregatedPSP = bool(int(sys.argv[6]))
    AnalyticShortestPath  = bool(int(sys.argv[7]))
    ShortestPathPrecuts = bool(int(sys.argv[8]))

if len(sys.argv) > 9:
    AnalyticCuts2 = bool(int(sys.argv[9]))
    LinkFiltering2 = bool(int(sys.argv[10]))


Stats = collections.defaultdict(int)
Stats['Alpha'] = alpha
Stats['Iterations'] = 0
Stats['LazyCuts'] = 0
Stats['Hubs'] = NHubs
Stats['AnalyticCuts'] = DisaggregatedAnalytic
Stats['PSP'] = DisaggregatedPSP
Stats['AnalyticSP'] = AnalyticShortestPath
Stats['PreCuts'] = ShortestPathPrecuts
Stats['Day'] = day
Stats['Heuristic'] = Heuristic
Stats['AnalyticCuts2'] = AnalyticCuts2
Stats['LinkFiltering2'] = LinkFiltering2
alpha = alpha/3600 #convert $/hr to $/sec
"""
Initialize sets and import data
"""

# hubs may be 10 or 20
if NHubs == 10:
    H = ImportHubs('10hubs.dat')
elif NHubs == 20:
    H = ImportHubs('20hubs-HOMEMADE.dat')
elif NHubs == 30:
    H = ImportHubs('30hubs-HOMEMADE.dat')
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

SHORTEST, LEGS, removed = TripsAndLegs(T,H,Tau,Gamma,Beta,D)


print('bundles of legs:',len(LEGS))

print(len(D),len(removed),len(D)+len(removed))


SHORTESTFWD = {}
for (Or,De),dr in D.items():
    SHORTESTFWD[Or,De] = DijkstraFwd(dr, Or, De, 
               {(h,l): 1 for (h,l) in Beta}
               , Tau,Gamma,H)
    
"""
Link Filtering
"""

if LinkFiltering:
    Links(D, H, Tau, Gamma)

if LinkFiltering2:
    usedBuses, unusedBuses = Links2(D, SHORTEST, SHORTESTFWD, Gamma, Tau)
    num = 0
    for key,List in unusedBuses.items():
        num += len(List)
    print(num)

"""
Set up BMP variables and objective
"""
start = time.time()
if DisaggregatedAnalytic or DisaggregatedPSP:
    BMP = Model("Benders Master Problem")
    BMP.setParam('MIPGap', 0.00)

    Z = {(h,l): BMP.addVar(vtype=GRB.BINARY)  for (h,l) in Gamma}
    Theta = {r: BMP.addVar() for r in D}
    Stats['SubProblems'] = len(D)
         
    BMP.setObjective(quicksum(Beta[h,l] * Z[h,l] for (h,l) in Beta) + 
                     quicksum(T[r]*Theta[r] for r in D)
                     , GRB.MINIMIZE)

    WeakConnectivity = {h:
        BMP.addConstr(quicksum(Z[h,l] for l in H if l != h) ==
                      quicksum(Z[l,h] for l in H if l != h)) for h in H}

  
                               
"""
Analytic Callback function
"""         
def  DisaggregatedAnalyticCallback(model, where):
    if where == GRB.Callback.MIPSOL:
        Stats['Iterations'] += 1
        ZBar = {k:v for (k,v) in zip(Z.keys(),model.cbGetSolution(Z.values()))}
        ThetaD = {k:v for (k,v) in zip(Theta.keys(),
                  model.cbGetSolution(Theta.values()))}
  
        for (Or,De),dr in D.items():
            sToDe = SHORTEST[Or,De]
            sFromOr = SHORTESTFWD[Or,De]
            toDe, xxxx = DijkstraRev(dr, Or, De, ZBar,Tau,Gamma,H)
            if AnalyticCuts2:
                fromOr = DijkstraFwd(dr, Or, De, ZBar,Tau,Gamma,H)
            SPath = toDe[Or]

            if SPath > ThetaD[Or,De] + EPS:
                V = {}
                if AnalyticCuts2:
                    before = {}
                    after = {}
                    for h in H:
                        if sFromOr[h] < fromOr[h]:
                            before[h] = 1
                        else :
                            before[h] = 0
                        if sToDe[h] < toDe[h]:    
                            after[h] = 1
                        else: 
                            after[h] = 0

                    for (h,l) in Beta:
                        if ZBar[h,l] > 0.5:
                            V[h,l] = 0
                        else:
                            step1 = SPath - fromOr[h] - toDe[l] - Gamma[h,l]
                            if before[h] == 1 or after[l] == 1:
                                step2 = (SPath - sFromOr[h] - sToDe[l]- Gamma[h,l])/2
                                V[h,l] = max(0,step1,step2)
                            else:
                                V[h,l] = max(0,step1)
                                   
                if not AnalyticCuts2:
                    U = {}
                    for h in H:
                        U[h] = max(0,min(SPath-sFromOr[h],toDe[h]))
                    for (h,l) in Beta:
                        if ZBar[h,l] > 0.5:
                            V[h,l] = 0
                        if ZBar[h,l] <= 0.5:
                            V[h,l] = max(0,U[h]- U[l] - Gamma[h,l])
                            
                model.cbLazy(Theta[Or,De]>=SPath - 
                             quicksum(Z[h,l] * V[h,l]
                             for (h,l) in Beta))

                Stats['LazyCuts'] += 1

                
"""
DUAL DISAGGREGATED BENDERS SUB PROBLEM
"""      
if not AnalyticShortestPath:
    BDSP = Model("Benders Dual of Sub Problem")
    BDSP.setParam('OutputFlag',0)
    
    UOr = BDSP.addVar()
    UDe = BDSP.addVar()
    U = {h: BDSP.addVar()  for h in H}
    V = {(h,l): BDSP.addVar() for (h,l) in Beta}
    
    ZBar = {k: 1 for k in Z}
    
    # RHS of contrainsts depend on trip Or/De
    OD = BDSP.addConstr(UOr - UDe <= 100)
    OH = {h: BDSP.addConstr(UOr - U[h] <= 100) for h in H}
    LD = {h: BDSP.addConstr(U[h] - UDe <= 100) for h in H}
    
    HL = {(h,l): BDSP.addConstr(U[h] - U[l] - V[h,l] <= Gamma[h,l])
            for (h,l) in Beta}

"""
PARETO SUB-PROBLEM
"""
if DisaggregatedPSP :
    PSP = Model("Pareto Sub Problem")
    PSP.setParam('OutputFlag',0)

    UOr2 = PSP.addVar()
    UDe2 = PSP.addVar()
    U2 = {h: PSP.addVar()  for h in H}
    V2 = {(h,l): PSP.addVar() for (h,l) in Beta}

    ZCore = {k: 0.1 for k in Z}

    # RHS of contrainsts depend on trip Or/De
    OD2 = PSP.addConstr(UOr2 - UDe2 <= 100)
    OH2 = {h: PSP.addConstr(UOr2 - U2[h] <= 100) for h in H}
    LD2 = {h: PSP.addConstr(U2[h] - UDe2 <= 100) for h in H}

    HL2 = {(h,l): PSP.addConstr(U2[h] - U2[l] - V2[h,l] <= Gamma[h,l])
            for (h,l) in Beta}

    global Equality
    Equality = PSP.addConstr(UOr2-UDe2 - 
            quicksum(1 * V2[h,l] for (h,l) in Beta) == 1 )

    PSP.update()
    
    
def  DisaggregatedPSPCallback(model, where):
    if where == GRB.Callback.MIPSOL:
        Stats['Iterations'] += 1
        ZBar = {k:v for (k,v) in zip(Z.keys(),model.cbGetSolution(Z.values()))}
        ThetaD = {k:v for (k,v) in zip(Theta.keys(),
                  model.cbGetSolution(Theta.values()))}
                     
        global AnalyticShortestPath
        if not AnalyticShortestPath :
            BDSP.setObjective(UOr - UDe -
                    quicksum(ZBar[h,l]*V[h,l] for (h,l) in V),
                    GRB.MAXIMIZE)
            
        for (h,l) in Beta:
            ZCore[h,l] = (ZCore[h,l]+ ZBar[h,l])/2
             
        PSP.setObjective(UOr2 - UDe2 -
            quicksum(ZCore[h,l]*V2[h,l] for (h,l) in V2),
            GRB.MAXIMIZE)
            
        global Equality
        PSP.remove(Equality)
        Equality = PSP.addConstr(UOr2-UDe2- 
            quicksum(ZBar[h,l] * V2[h,l] for (h,l) in Beta) == 1 )
            
        for (Or,De),dr in D.items() :
            if AnalyticShortestPath :
                SPath = DijkstraShortestPath(dr, Or, De, ZBar,Tau,Gamma,H)
#                SPath = AStar(dr, Or, De, ZBar,Tau,Gamma,H,SHORTEST[Or,De])
                    
            if not AnalyticShortestPath :
                OD.rhs = Tau[Or,De]
                for h in H:
                    OH[h].rhs = Tau[Or,h]
                    LD[h].rhs = Tau[h,De] 
                        
                BDSP.optimize()
                SPath = BDSP.objVal
                      
            if SPath > ThetaD[Or,De]:
                Equality.rhs = SPath
                OD2.rhs = Tau[Or,De]
                for h in H:
                    OH2[h].rhs = Tau[Or,h]
                    LD2[h].rhs = Tau[h,De] 
                                  
                PSP.optimize()

                if PSP.Status == 4:
                    PSP.setParam('OutputFlag',1)
                    PSP.setParam('DualReductions',0)
                    PSP.optimize()
                    print(PSP.Status)
                    PSP.setParam('DualReductions',1)
                    PSP.setParam('OutputFlag',0)

                Stats['LazyCuts'] += 1
                model.cbLazy(Theta[Or,De]>=
                            UOr2.x - UDe2.x -
                            quicksum(Z[h,l] * V2[h,l].x for (h,l) in Beta))


if ShortestPathPrecuts: # Do precuts
    for (Or,De),dr in D.items():
        BMP.addConstr(Theta[Or,De] >= SHORTEST[Or,De][Or])   

loadtime = time.time() - start

Stats['LoadTime'] = loadtime
    
if DisaggregatedAnalytic:
    BMP.setParam('LazyConstraints',1)
    BMP.setParam('Threads',4)
    BMP.setParam('TimeLimit',1200.0)
    BMP.optimize(DisaggregatedAnalyticCallback)
    Stats['ObjVal'] = BMP.objVal
    
if DisaggregatedPSP:
    BMP.setParam('LazyConstraints',1)
    BMP.setParam('Threads',4)
    BMP.setParam('TimeLimit',1200.0)
    BMP.optimize(DisaggregatedPSPCallback)
    Stats['ObjVal'] = BMP.objVal

Stats['RunTime'] = time.time() - start
Stats['Trips'] = len(D)
#for z in Z:
#    if Z[z].x > 0.9:
#        print(Z[z].x, z)

if len(sys.argv) > 1:
    with open(str(sys.argv[1]), "a") as f:
        f.write(sys.argv[0]+','+str(datetime.datetime.now())+','+
                str(sorted([(k,v) for (k,v) in zip(Stats.keys(),Stats.values())]))+'\n')  
