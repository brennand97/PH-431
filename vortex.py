import numpy as np
import matplotlib.pyplot as plt
from scipy import special as sp

mu = 1
e0 = 1
c = 1
w = 1
z = 1
l = 4
p = 2
A0 = 1

def rho(x,y):
    return np.sqrt(x**2+y**2) 

def u(x, y, l, p):
    r = rho(x,y)
    return (-1)**p * (np.sqrt(2) * r/w)**l * \
        sp.genlaguerre(p, l)(2 * r**2 / w**2) * \
        np.exp(- r**2 / w**2)

def I(x, y, l, p): 
    return 0.5 / (mu * c) * A0**2 * ( u(x, y, l, p) )**2

x = np.linspace(-4,4,500)
y = x
X,Y = np.meshgrid(x,y) 
plt.figure()
plt.pcolor(X,Y, I(X,Y,l,p), cmap='gray')

plt.title(r"Optical vortex intensity cross section z={}".format(z))

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

plt.show()