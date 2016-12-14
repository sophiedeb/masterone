#-------------------------------------------------------------------------------
#                import necessary packages
#-------------------------------------------------------------------------------

import multiprocessing
import numpy as np
from collections import Counter
import random

num_cores = multiprocessing.cpu_count()        # 4 on quad-core pc, 12 on Helios

#-------------------------------------------------------------------------------
#                define the necessary parameters
#-------------------------------------------------------------------------------

S = 25                  # number of species in the metacommunity pool
r = 20                  # number of species in each local community
LCs = 500               # number of local communities (LCs)
q = 50000               # number of lattice sites
connectance = 0.9       # fraction of nonzero elements in the interaction matrix

SIS = [10,15,20]        # list of strongly interacting species (SISs)
SISfactor = 15          # relative strength of SISs
minmigration = 0.1
maxmigration = 0.1
minextinction = 0.04
maxextinction = 0.04

# convergence criterion: relaxed for the empty sites
epsilon = [0.01*q for _ in range(S+1)]
epsilon[0] = 0.2*q

#-------------------------------------------------------------------------------
#                help functions
#-------------------------------------------------------------------------------

# check if convergence has been reached ----------------------------------------

# the main idea is that the occurence of each individual species shouldn't have
# changed more than a predefined constant epsilon during the last 200 time steps
def Convergence(statevectors,timestep):
    if timestep <= 211 or timestep%10 != 0:
        return False
    else:
        check = True
        for t in range(0,200/10):
            diff = np.subtract(statevectors[-1],statevectors[-(t+2)])
            check = (diff < np.array(epsilon)).all()
            if check == False: break
        return check

# choose a random site, the species on it shouldn't have changed in this timestep
def randomSite(lattice, templattice):
    k2 = 0
    while True:
        k2 = random.randrange(0,q,1)
        if lattice[k2] == templattice[k2]: break
    j = lattice[k2]
    return k2, j

#-------------------------------------------------------------------------------
#                initialize stuff
#-------------------------------------------------------------------------------

# interaction matrix -----------------------------------------------------------

omega = np.zeros((S+1, S+1))

# fill up the "background" interaction
for row in range(1,S+1):
    for column in range(S+1):
        if row != column: omega[row][column] = random.uniform(-1,1)/SISfactor

# SISs are defined by high probabilities in their column
for species in SIS:
    omega[:,species] = [random.uniform(-1,1) for _ in range(S+1)]
    omega[species][species] = 0

# the growth rates are uniformly distributed from [0,1)
omega[:,0] = [random.uniform(0,1) for _ in range(S+1)]
omega[0,:] = [0 for _ in range(S+1)]

# apply the "connectance" (fraction of nonzero elements in omega)
for matrixelement in range(int((1-connectance)*(S)**2)):
    while True:
        (row, column) = random.randrange(1,S+1,1), random.randrange(1,S+1,1)
        if omega[row][column] != 0:
            omega[row][column] = 0
            break

# immigration rates ------------------------------------------------------------

mu = [random.uniform(minmigration,maxmigration) for _ in range(S+1)]
mu[0] = 0

# extinction rates -------------------------------------------------------------

e = [random.uniform(minextinction,maxextinction) for _ in range(S+1)]
e[0] = 0

# save data for later use or reference -----------------------------------------

np.savetxt("omega.txt",omega)
np.savetxt("mu.txt",mu)
np.savetxt("e.txt",e)
np.savetxt("SIS.txt",SIS)

#-------------------------------------------------------------------------------
#                the actual loop function
#-------------------------------------------------------------------------------

def Loop(s):
    # choose r species from the pool of S:
    specieslist = []
    while len(specieslist) < r:
        trial = random.randrange(1,S+1,1)
        if trial not in specieslist:
            specieslist.append(trial)

    # initialize the lattice
    lattice = [0 for _ in range(q)]
    for i in specieslist:
        for k in range(random.randrange(1,q/r,1)):
            lattice[random.randrange(0,q,1)] = i

    timestep = 0
    statevectors = []

    # loop in time until we reach convergence
    while not Convergence(statevectors,timestep):
        for k in range(q):                       # loop over the lattice sites
            A = lattice[k]                       # species on lattice site k
            B = random.choice(specieslist)
            if A == 0 and random.uniform(0,1) <= mu[B]:
                lattice[k] = B                   # migration
            elif A != 0 and random.uniform(0,1) <= e[A]:
                lattice[k] = 0                   # extinction

        # a species should only interact once every timestep, otherwise it can
        # do weird stuff like dying more than once etc. That's why we introduce
        # a constant reference lattice that doesn't change during 1 timestep
        templattice = lattice
        for k1 in range(q):
            if templattice[k1] != lattice[k1]: continue
            A = lattice[k1]
            k2, B = randomSite(lattice, templattice)
            probability = omega[B][A] + np.absolute(omega[A][B])

            if omega[A][B] < 0 and random.uniform(0,1) <= probability:
                lattice[k1] = B

            elif omega[A][B] > 0 and omega[B][A] >= 0:
                k3, C = randomSite(lattice, templattice)
                if probability > omega[C][A] and random.uniform(0,1) <= probability:
                    lattice[k3] = A
                # no distinction was made if C == 0 or not, because by
                # construction we have omega[C][A] = 0 for all A

        # determine the state vector at a given timestep and add to the statevectors
        # NOTE: counting is very time consuming, don't do it every timestep
        # NOTE: the count() function is O(n^2), Counter() is O(n)
        if timestep%10 == 0:
            c = Counter(lattice)
            x = np.zeros(S+1)
            for i in range(0,S+1):
                x[i] = c[i]
            statevectors.append(x)

        timestep += 1

        if timestep == 6000:
            print "Finished " + str(s) + "/" + str(LCs) + " without reaching convergence"
            statevectors = [[np.nan for _ in range(S+1)]]
            # when the program is finished delete the nan's manually with CTRL+F
            break

    # keep track of progress
    print "Finished " + str(s) + "/" + str(LCs) + " after " + str(timestep-1) + " iterations"

    # the equilibrium statevector is the final statevector
    return statevectors[-1]

#-------------------------------------------------------------------------------
#            execute Loop() for 500 LCs, using parallel programming
#-------------------------------------------------------------------------------

#equilibria = []
#for s in range(1,LCs+1):
#    equilibria.append(Loop())

pool = multiprocessing.Pool(num_cores)
equilibria = pool.map(Loop, range(1, LCs+1, 1))

#-------------------------------------------------------------------------------
#           save the data, to be used in NMDS.R
#-------------------------------------------------------------------------------

np.savetxt("equilibria.txt", equilibria)

print "Done!"
