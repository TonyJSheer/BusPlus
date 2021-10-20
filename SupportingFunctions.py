



def ImportHubs(filename):
    H = []
    File1 = open(filename,'r')
    for line in File1:
        H.append(int(line[0:len(line)-1]))
    File1.close()
    return H


def BuildCostsOriginal(filename,c,b,S,n,alpha,H):
    Tau = {}
    Gamma = {}
    Beta = {}
    N = [4026]
    
    File2 = open(filename,'r')
    for num, line in enumerate(File2):
        if num != 0:
            i,j,Dist,Time = line.split(',')
            i = int(i)
            j = int(j)
            Tau[i,j] = (1 - alpha)*c* (float(Dist)) + alpha * float(Time)
            if i in H and j in H and i != j:
                Gamma[i,j] = alpha * (float(Time) + S)
                Beta[i,j] = (1-alpha) * b * n * float(Dist)
            if i == 4026:
                N.append(j)
    for i in N:
        Tau[i,i] = 0
    File2.close()
    return N,Tau,Gamma,Beta

def BuildCosts(filename,c,b,S,n,alpha,H):
    Tau = {}
    Gamma = {}
    Beta = {}
    N = [4026]
    
    File2 = open(filename,'r')
    for num, line in enumerate(File2):
        if num != 0:
            i,j,Dist,Time = line.split(',')
            i = int(i)
            j = int(j)
            Tau[i,j] = c* (float(Dist)) + alpha * float(Time)
            if i in H and j in H and i != j:
                Gamma[i,j] = alpha * (float(Time) + S)
                Beta[i,j] = b * n * float(Dist)
            if i == 4026:
                N.append(j)
    for i in N:
        Tau[i,i] = 0
    File2.close()
    return N,Tau,Gamma,Beta


def ImportTrips(filename,Tau):
    T = {}
    File3 = open(filename,'r')
    for num, line in enumerate(File3):
        if num != 0:
            Time, Or, De = line.split(',')
            if Or != 'Unknown' and De[0:len(De)-1] != 'Unknown':
                Or = int(Or)
                De = int(De)
                if Or!=De and (Or,De) in Tau:
                    try:
                        T[Or,De] +=1 #if its already there, the trip happens again
                    except KeyError:
                        T[Or,De]=1
    File3.close()
    return T

def TripsAndLegs(T,H,Tau,Gamma,Beta,D):
    SHORTEST = {}
    LEGS = {}
    removed=[]
    for Or,De in T:
        toDe, prevDe = DijkstraRev(D[Or,De], Or, De,
                             {(h,l): 1 for (h,l) in Beta},
                             Tau,Gamma,H)
        SHORTEST[Or,De] = toDe
        if toDe[Or] > Tau[Or,De]:
            crash
        if toDe[Or] == Tau[Or,De]:
            del D[Or,De]
            removed.append((Or,De))
        if toDe[Or] < Tau[Or,De]:
            # Check for degenerate trips(no bus arcs)
            currentNode = Or
            while currentNode != De:
                parentNode, arcType = prevDe[currentNode]
                if arcType == 'b' :
                    busLegStart = currentNode
                    busLegEnd = parentNode
                    try :
                        LEGS[busLegStart,busLegEnd].append((Or,De))
                    except KeyError:
                        LEGS[busLegStart,busLegEnd] = [(Or,De)]
                    Gamma[busLegStart, busLegEnd] # Will crash if not valid bus arc
                    break
                currentNode = parentNode
            if arcType == 't':
                del D[Or,De]
                removed.append((Or,De))
    return SHORTEST, LEGS, removed


def Links(D, H, Tau, Gamma):
    for (Or,De),dr in D.items():
        for h in H:
            if (h != Or and h != De):
                arcUsedOr = False
                arcUsedDe = False
                for l in H :
                    if l != h:
                        if Tau[Or,h] + Gamma[h,l] + Tau[l,De] <= Tau[Or,De]:
                            arcUsedOr = True
                        if Tau[Or,l] + Gamma[l,h] + Tau[h,De] <= Tau[Or,De]:
                            arcUsedDe = True
                if Tau[Or,h] + Tau[h,De] <= Tau[Or,De]:
                    arcUsedOr = True
                    arcUsedDe = True
                if not arcUsedOr:
                    dr.remove((Or,h))
                if not arcUsedDe:
                    dr.remove((h,De))
    return                      

def Links2(D, SHORTEST, SHORTESTFWD, Gamma, Tau):
    unusedBuses = {}
    usedBuses = {}
    for (Or,De) in D:
        unusedBuses[Or,De] = []
        usedBuses[Or,De] = []
        sFromOr, sToDe = SHORTESTFWD[Or,De], SHORTEST[Or,De]
        direct = Tau[Or,De]
        for (h,l),arcCost in Gamma.items():
            if sFromOr[h] + sToDe[l] + arcCost >= direct:
                unusedBuses[Or,De].append((h,l))
            else: 
                usedBuses[Or,De].append((h,l))
    return usedBuses, unusedBuses
    
def DijkstraFwd(dr, Or, De, ZBar, Tau,Gamma,H):
    fromOr = {}
    fromOr[Or] = 0
    ToExpand = []
    PrevOr = {}
#    # Get Or successors, add to PQ
    for (i,j) in dr:
        if i == Or:
            ToExpand.append((Tau[Or,j],j,Or))
    if Or in H:
        for i,(dist,dest,parent) in enumerate(ToExpand):
            if dest in H and Gamma[Or,dest] < dist and ZBar[Or,dest] > 0.9:
                ToExpand.pop(i)
                ToExpand.insert(i,(Gamma[Or,dest],dest,Or))

    while fromOr.__len__() < (dr.__len__()+3)/2 : 
        ToExpand.sort()
        (d,h,p) = ToExpand.pop(0)
        fromOr[h] = d
        PrevOr[h] = p
        if h != De:
            for i,(dist,dest,parent) in enumerate(ToExpand):
                if dest not in fromOr and dest == De:
                    if Tau[h,De] + d < dist:
                        ToExpand.pop(i)
                        ToExpand.insert(i,(Tau[h,De] + d,dest,h))
                    if De in H and ZBar[h,De] > 0.9 and \
                    Gamma[h,De] + d < dist and Gamma[h,De] < Tau[h,De]:
                        ToExpand.pop(i)
                        ToExpand.insert(i,(Gamma[h,dest] + d,dest,h))
                if dest not in fromOr and dest != De:
                    if Gamma[h,dest] + d < dist and ZBar[h,dest] > 0.9:
                        ToExpand.pop(i)
                        ToExpand.insert(i,(Gamma[h,dest] + d,dest,h))     
    return fromOr


def DijkstraRev(dr, Or, De, ZBar, Tau,Gamma,H):
    """
    Returns a shortest path tree from each node to destination(toDe).
    Also returns the previous node and arc type
    for each node that isn't the Destination
    """
    # Reverse
    toDe = {}
    toDe[De] = 0
    ToExpand = []
    prevDe = {}
    # Get De successors, add to PQ2
    for (i,j) in dr:
        if j == De:
            ToExpand.append((Tau[i,j],i,j,'t'))
    if De in H:
        for i,(dist,h,parent,TYPE) in enumerate(ToExpand):
            if h in H and Gamma[h,De] < dist and ZBar[h,De] > 0.9:
                ToExpand.pop(i)
                ToExpand.insert(i,(Gamma[h,De],h,De,'b'))

    while toDe.__len__() < (dr.__len__()+3)/2 : 
        ToExpand.sort()
        (d,h,p,t) = ToExpand.pop(0)
        toDe[h] = d
        prevDe[h] = (p,t)
        if h != Or:
            for i,(dist,dest,parent,TYPE) in enumerate(ToExpand):
                if dest not in toDe and dest == Or:
                    if Tau[Or,h] + d < dist:
                        ToExpand.pop(i)
                        ToExpand.insert(i,(Tau[Or,h] + d,Or,h,'t'))
                    if Or in H and ZBar[Or,h] > 0.9 and \
                    Gamma[Or,h] + d < dist and Gamma[Or,h] < Tau[Or,h]:
                        ToExpand.pop(i)
                        ToExpand.insert(i,(Gamma[Or,h] + d,Or,h,'b'))
                if dest not in toDe and dest != Or:
                    if Gamma[dest,h] + d < dist and ZBar[dest,h] > 0.9:
                        ToExpand.pop(i)
                        ToExpand.insert(i,(Gamma[dest,h] + d,dest,h,'b'))   

    return toDe, prevDe

def DijkstraShortestPath(dr, Or, De, ZBar, Tau,Gamma,H):
    # Reverse
    toDe = {}
    toDe[De] = 0
    ToExpand = []
    PrevDe = {}
    # Get De successors, add to PQ2
    for (i,j) in dr:
        if j == De:
            ToExpand.append((Tau[i,j],i,j))
    if De in H:
        for i,(dist,h,parent) in enumerate(ToExpand):
            if h in H and Gamma[h,De] < dist and ZBar[h,De] > 0.9:
                ToExpand.pop(i)
                ToExpand.insert(i,(Gamma[h,De],h,De))

    while toDe.__len__() < (dr.__len__()+3)/2 : 
        ToExpand.sort()
        (d,h,p) = ToExpand.pop(0)
        toDe[h] = d
        PrevDe[h] = p
        if h == Or:
            return d
        if h != Or:
            for i,(dist,dest,parent) in enumerate(ToExpand):
                if dest not in toDe and dest == Or:
                    if Tau[Or,h] + d < dist:
                        ToExpand.pop(i)
                        ToExpand.insert(i,(Tau[Or,h] + d,Or,h))
                    if Or in H and ZBar[Or,h] > 0.9 and \
                    Gamma[Or,h] + d < dist and Gamma[Or,h] < Tau[Or,h]:
                        ToExpand.pop(i)
                        ToExpand.insert(i,(Gamma[Or,h] + d,Or,h))
                if dest not in toDe and dest != Or:
                    if Gamma[dest,h] + d < dist and ZBar[dest,h] > 0.9:
                        ToExpand.pop(i)
                        ToExpand.insert(i,(Gamma[dest,h] + d,dest,h))   

    return toDe

def AStar(dr, Or, De, ZBar, Tau,Gamma,H, Heur):
    fromOr = {}
    fromOr[Or] = 0
    ToExpand = []
    PrevOr = {}
#    # Get Or successors, add to PQ
    for (i,j) in dr:
        if i == Or:
            cost = Tau[Or,j] 
            ToExpand.append((cost + Heur[j],j,Or, cost))
    if Or in H:
        for i,(distH,dest,parent,dist) in enumerate(ToExpand):
            if dest in H and Gamma[Or,dest] < dist and ZBar[Or,dest] > 0.9:
                ToExpand.pop(i)
                ToExpand.insert(i,(Gamma[Or,dest] + Heur[dest],
                                   dest,Or, Gamma[Or,dest]))

    while fromOr.__len__() < (dr.__len__()+3)/2 : 
        ToExpand.sort()
        (d,h,p,c) = ToExpand.pop(0)
        fromOr[h] = c
        PrevOr[h] = p
        if h == De:
            return c
        if h != De:
            for i,(distH,dest,parent,dist) in enumerate(ToExpand):
                if dest not in fromOr and dest == De:
                    arcCost = Tau[h,De]
                    if arcCost + c < dist:
                        ToExpand.pop(i)
                        ToExpand.insert(i,(arcCost + c + Heur[De],
                                           dest,h, arcCost + c))
                        
                    if De in H and ZBar[h,De] > 0.9 and \
                    Gamma[h,De] + c < dist and Gamma[h,De] < arcCost:
                        ToExpand.pop(i)
                        ToExpand.insert(i,(Gamma[h,De] + c + Heur[De],
                                           dest,h,Gamma[h,De] + c))
                        
                if dest not in fromOr and dest != De:
                    if Gamma[h,dest] + c < dist and ZBar[h,dest] > 0.9:
                        ToExpand.pop(i)
                        ToExpand.insert(i,(Gamma[h,dest] + c + Heur[dest],
                                           dest,h, Gamma[h,dest] + c))     
    return fromOr