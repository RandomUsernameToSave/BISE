from BISE import Species,Reaction,Solver


Te = 300
Tg = 300
P = 1e5

kb = 1.380649e-23
n = P / kb/ Tg

Xe =      Species('Xe',n)
e =       Species('e', 1e4)
Xep =     Species('Xe^+', 1e4)
Xe2p =    Species('Xe2^+', 1.)
Xe3p =    Species('Xe3^+', 1.)
Xess =    Species('Xe**', 1.)
Xe2ss =   Species('Xe2**',1.)
Xes =     Species('Xe*',1.)
Xe2s_1 =  Species('Xe2*(1)', 1.)
Xe2s_3 =  Species('Xe2*(3)', 1.)

R1 = Reaction([Xep, Xe, Xe] , [Xe2p, Xe], 2*1e-31 * (300/Tg)**0.5)
R2 = Reaction([Xe2p, Xe , Xe] , [Xe3p,Xe], 6 * 1e-32 * (300/Tg)**0.5)
R3 = Reaction([Xe2p,e], [Xess,Xe], 2.3*1e-7 * Te**(-0.7))
R4 = Reaction([Xe3p, e], [Xess,Xe,Xe] , 1e-5 * Te**(-0.5))
R5 = Reaction([Xess,Xe,Xe], [Xe2ss,Xe] , 5e-31*(300/Tg)**0.5)
R6 = Reaction([Xe2ss], [Xes,Xe],1e8)
R7 = Reaction([Xe2ss,Xe] , [Xes,Xe,Xe], 1e-11 * (Tg/300)**0.5)
R8 = Reaction([Xes, Xe,Xe] , [Xe2s_3,Xe] , 4.4e-32 * (300/Tg)**0.5)
R9 = Reaction([Xes,Xe,Xe], [Xe2s_1,Xe], 2e-32 * (300/Tg)**0.5)
R10a = Reaction([e,Xe2s_3],[e,Xe2s_1], 1.8e-8)
R10b = Reaction([e,Xe2s_1],[e,Xe2s_3], 4.9e-8)
R11a = Reaction([Xe,Xe2s_3],[Xe,Xe2s_1], 4.6e-15 * (Tg/300)**0.5)
R11b = Reaction([Xe,Xe2s_1],[Xe,Xe2s_3], 1.2e-13 * (Tg/300)**0.5)
R12 = Reaction([Xe2s_1], [Xe,Xe] , 2.1e8)
R13 = Reaction([Xe2s_3], [Xe,Xe] , 1e7)
R16a = Reaction([Xe2s_1,Xe2s_1] , [Xe2p,Xe,Xe,e], 8e-11)
R16b = Reaction([Xe2s_3,Xe2s_3] , [Xe2p,Xe,Xe,e], 8e-11)
R17 = Reaction([Xes,Xes], [Xep,Xe,e],8e-11)
R18a = Reaction([Xe2s_1,e], [Xe2p,e,e], 5e-9)
R18b = Reaction([Xe2s_3,e], [Xe2p,e,e], 5e-9)
R19 = Reaction([Xes,e], [Xep,e,e] , 2.7e-9 )
R20a = Reaction([Xe2s_1,e], [Xe,Xe,e],4e-9)
R20b = Reaction([Xe2s_3,e], [Xe,Xe,e],4e-9)
R21 = Reaction([Xes,e], [Xe,e] , 3e-9 )
R22a = Reaction([Xe2s_1,e], [Xes,Xe,e], 2e-7)
R22b = Reaction([Xe2s_3,e], [Xes,Xe,e], 2e-7)
R23a = Reaction([Xe2s_1,e], [Xe2ss,e] , 3e-7)
R23b = Reaction([Xe2s_3,e], [Xe2ss,e] , 3e-7)
R24a = Reaction([Xe2ss,e], [Xe2s_3,e],6e-7)
R24b = Reaction([Xe2ss,e], [Xe2s_1,e],2e-7)
R25 = Reaction([Xess], [Xes], 1e-6)
R26 = Reaction([Xess,e], [Xep,e,e], 2e-8)
R27 = Reaction([Xes,e], [Xess,e], 3e-7)
R28 = Reaction([Xess,e], [Xes,e], 8e-7)
R30 = Reaction([Xess,Xess], [Xep,Xe,e], 1e-10)
R31 = Reaction ([Xe2ss],[Xe2p , Xe,Xe,e] ,1e-10)
R32 = Reaction([Xe2p,e], [Xep,Xe,e],1e-7)
R33 = Reaction([Xe2ss,e], [Xe2p,e,e],6e-8)

solver = Solver()

solver.add_species([Xe,Xe2p,Xe2s_1,Xe2s_3,Xe2ss,Xe3p,Xes,Xess,e,Xep])
solver.add_reactions([R1,R2,R3,R4,R5,R6,R7,R8,R9,R10a,R10b,R11a,R11b,R12,R13,R16a,R16b,R17,R18a,R18b,R19,R20a,R20b
                      ,R21,R22a,R22b,R23a,R23b,R24a,R25,R26,R27,R28,R30,R31,R32,R33])
solver.summary()
solver.init()
sol = solver.solve([0, 3e-5])
import matplotlib.pyplot as plt
#plt.plot(sol.y[0])
plt.plot(sol.y[1])
plt.plot(sol.y[2])
plt.show()