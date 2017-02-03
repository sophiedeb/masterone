import numpy as np
import random
from scipy.integrate import odeint

S = 50
s = 40
LCs = 500
SIS = [20]
SISfactor = [170]

r = [random.uniform(0,1) for _ in range(S)]

omega = np.array([[random.uniform(-0.03,0.03) for _ in range(S)] for _ in range(S)])

for row in range(0,S):
    omega[row][row] = -1

index = 0
for species in SIS:
    omega[:,species] *= SISfactor[index]
    omega[species][species] = -1
    index += 1

t = np.linspace(0,100,2)

equilibria = []

def loop(l):
    specieslist = []
    for i in range(s):
        while True:
            trial = random.randrange(S)
            if trial not in specieslist:
                specieslist.append(trial)
                break
    
    def f(x,t):
        dxdt = np.zeros(S)
        for i in specieslist:
            dxdt[i] = r[i]*x[i]
            for j in specieslist:
                dxdt[i] += x[i]*omega[i][j]*x[j]
            if dxdt[i] < 0: dxdt[i] = 0
        return dxdt

    x0 = np.array([0. for _ in range(S)])
    for i in specieslist:
        x0[i] = random.uniform(0,1)

    sol = odeint(f,x0,t)
    x = np.array(sol[-1])

    if (np.linalg.norm(f(x,0)) < 0.01):
        print "CONVERGENCE REACHED!"
        for i in range(S):
            if x[i] < 0: x[i] = 0
        equilibria.append(x)
    else:
        print "CONVERGENCE NOT REACHED!"

for l in range(LCs):
    loop(l)

np.savetxt("equilibria20.txt", equilibria)
np.savetxt("SIS.txt", SIS)
