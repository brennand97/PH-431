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
<<<<<<< Updated upstream
=======
    """
    Evaluates the amplitude of a helical beam, at coordinates x, y relative to the
	beam's axis using equation 1 from Sundbeck; a special case of "Laguerre-Gaussian
	eigenmodes of the Helmholtz equation."

    Parameters:
    x, y: float
    Cartesian coordinates in a 2D plane, where the optical axis normal to the plane.

    l: int
    Azimuthal index of the generalized Laguerre polynomial. Also called the helicity,
	or topological charge of the beam.

    p: int
    Radial index of the generalized Laguerre polynomial.

    Returns:
    u: float
    Real-valued amplitude of a helical beam at coordinates x, y.

    References:
    Steven Sundbeck, Ilya Gruzberg, and David G. Grier, "Structure and scaling of
	helical modes of light," Opt. Lett. 30, 477-479 (2005)
    """

    # Helical beam has a radially symmetrical amplitude,
    # so the amplitude function is only dependent on the
    # distance from the origin to the x, y coordinates.
>>>>>>> Stashed changes
    r = rho(x,y)
    return (-1)**p * (np.sqrt(2) * r/w)**l * \
        sp.genlaguerre(p, l)(2 * r**2 / w**2) * \
        np.exp(- r**2 / w**2)

def I(x, y, l, p): 
<<<<<<< Updated upstream
=======
    """
    Evaluates the electromagnetic intensity of a helical beam at coordinates x, y
	relative to the beam's axis.

    Parameters:
    x, y: float
    Cartesian coordinates in a 2D plane, where the optical axis normal to the plane.

    l: int
    Azimuthal index of the generalized Laguerre polynomial. Also called the helicity,
	or topological charge of the beam.

    p: int
    Radial index of the generalized Laguerre polynomial.

    Returns:
    I: float
    Electromagnetic intensity of a helical beam at coordinates x, y.
    """
    
>>>>>>> Stashed changes
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