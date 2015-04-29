import truss
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import patches

modulus = 210.e3     #Pa
rho = 2700.           #kg/m**3
surface = .001       #m**2
yield_stress = 400 #Pa
from scipy import optimize

m = truss.core.Model()
A = m.add_node((0., 1.), label = "A")
B = m.add_node((.5, 1.), label = "B")
C = m.add_node((1., 1.), label = "C")
D = m.add_node((.5, 0.), label = "D")
E = m.add_node((1., 0.), label = "E")
F = m.add_node((2., 0.), label = "F")
G = m.add_node((2., 1.), label = "G")





A.block[1] = True
A.block[0] = True
C.block[0] = True
C.block[1] = True
D.block[1] = True
D.block[0] = True


AB = m.add_bar(A, B, modulus = modulus, density = rho, section = surface)
BC = m.add_bar(B, C, modulus = modulus, density = rho, section = surface)
BD = m.add_bar(B, D, modulus = modulus, density = rho, section = surface)

m.add_bar(C, E, modulus = modulus, density = rho, section = surface)
m.add_bar(E, F, modulus = modulus, density = rho, section = surface)
m.add_bar(F, G, modulus = modulus, density = rho, section = surface)

B.force  = np.array([-1., -1.])


#m.solve()

xlim, ylim = m.bbox(deformed = False)
fig = plt.figure(0)
plt.clf()
ax = fig.add_subplot(1,1,1)
ax.set_aspect("equal")
#ax.axis("off")
m.draw(ax, deformed = True, field = "stress", label = True, force_scale = .1, forces = True)
plt.xlim(xlim)
plt.ylim(ylim)
plt.grid()  


plt.show()
