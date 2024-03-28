
from BISE import Species,Reaction,Solver

Ar = Species('Ar', 2.5e19)
Arp = Species('Ar+', 1e0)
e = Species('e', 1e0)

R1 = Reaction([e,Ar] , [e,e, Arp], 1e-13,True)
R2 = Reaction([e,Arp,Ar] , [Ar, Ar],1e-25)

solver = Solver()
solver.add_species([Ar,Arp,e])
solver.add_reaction(R1)
solver.add_reaction(R2)
solver.init()
sol = solver.solve([0, 3e-2])
print(sol)
import matplotlib.pyplot as plt
#plt.plot(sol.y[0])
plt.plot(sol.y[1])
plt.plot(sol.y[2])
plt.show()


