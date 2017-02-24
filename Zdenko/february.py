"""
Implementation of the individual-based model discussed in February
@author: zheyvaer
"""

#-------------------------------------------------------------------------------
#                import necessary packages
#-------------------------------------------------------------------------------

import multiprocessing
import numpy as np
import random

num_cores = multiprocessing.cpu_count()        # 4 on quad-core pc, 12 on Helios

#-------------------------------------------------------------------------------
#                define the necessary parameters
#-------------------------------------------------------------------------------

S = 100                              # number of species in the metacommunity
s = 80                               # number of species in each local community
LCs = 500                            # number of local communities
sites = 5000                         # number of lattice sites
connectance = 1                      # the connectance

SIS = [20]
SISfactor_ = 200
maxint = 0.9
minmigration = 0.4
maxmigration = 0.4
minextinction = 0.3
maxextinction = 0.3

SISfactor = [1 for _ in range(S+1)]
for sis in SIS:
    SISfactor[sis] = SISfactor_

#-------------------------------------------------------------------------------
#                construct the interaction matrix
#-------------------------------------------------------------------------------

omega = [[random.uniform(-maxint,maxint) for _ in range(S+1)] for _ in range(S+1)]

for i in range(1,S+1):
    for j in range(1,S+1):
        if random.uniform(0,1) <= 1.-connectance:
            omega[i][j] = 0

for i in range(S+1):
    omega[i][0] = random.uniform(0,1)
    omega[0][i] = 0
    omega[i][i] = -1.

# immigration rates
mu = [random.uniform(minmigration,maxmigration) for _ in range(S+1)]
mu[0] = 0

# extinction rates
e = [random.uniform(minextinction,maxextinction) for _ in range(S+1)]
e[0] = 0

# save everything for later use and/or reference
np.savetxt("omega.txt",omega)
np.savetxt("mu.txt",mu)
np.savetxt("e.txt",e)
np.savetxt("SIS.txt",SIS)
np.savetxt("SISfactor.txt",SISfactor)

#-------------------------------------------------------------------------------
#           the loop describes the interaction in one community
#-------------------------------------------------------------------------------

def Loop(u):
    # choose s species from the pool of S:
    specieslist = []
    while len(specieslist) < s:
        trial = random.randrange(1,S+1,1)
        if trial not in specieslist:
            specieslist.append(trial)

    # initialize the lattice
    lattice = [0 for _ in range(sites)]
    occurenceprob = [random.uniform(0,1) for _ in range(S+1)]
    for site in range(sites):
        species = random.choice(specieslist)
        # we want to have at least 66% empty sites
        if random.uniform(0,1) <= 0.33 and random.uniform(0,1) <= occurenceprob[species]:
            lattice[site] = species

    # initialize the statevector
    x = [0 for _ in range(S+1)]
    for i in range(S+1):
        x[i] = lattice.count(i)

    # I don't have an actual convergence criterium due to the big fluctuations,
    # but 250 is always enough (most of the time even 50 is enough)
    for t in range(250):
        # initialize the so-called "interaction vector"
        y = [0 for _ in range(S+1)]
        for i in range(S+1):
            y[i] = x[i]*SISfactor[i]

        # the number of interactions
        q = sum(y)

        # temporary state vector
        xtemp = x

        # "lattice" from which we will draw random species A
        xlatt = []
        for species in range(S+1):
            for i in range(x[species]):
                xlatt.append(species)

        # "lattice" from which we will draw random species B
        ylatt = []
        for species in range(S+1):
            for i in range(y[species]):
                ylatt.append(species)

        for interaction in range(q):
            A = random.choice(xlatt)
            B = random.choice(ylatt)
            ext = True
            if omega[A][B] < 0 and random.uniform(0,1) <= abs(omega[A][B]):
                if xtemp[A] > 0:
                    xtemp[0] += 1
                    xtemp[A] -= 1
                    ext = False
            elif omega[A][B] > 0 and random.uniform(0,1) <= omega[A][B]:
                if xtemp[0] > 0:
                    xtemp[0] -= 1
                    xtemp[A] += 1

            if A == 0:
                C = random.choice(specieslist)
                if random.uniform(0,1) <= mu[C]:
                    xtemp[0] -= 1
                    xtemp[C] += 1
            else:
                if random.uniform(0,1) <= e[A] and xtemp[A] > 0 and ext:
                    xtemp[0] += 1
                    xtemp[A] -= 1

        x = xtemp

    print("Finished LC nr. " + str(u) + "/" + str(LCs))
    return x

#-------------------------------------------------------------------------------
#           save the data, to be used in NMDS.R
#-------------------------------------------------------------------------------

pool = multiprocessing.Pool(num_cores)
equilibria = pool.map(Loop, range(1, LCs+1, 1))

np.savetxt("equilibria.txt", equilibria)

print("Done!")
