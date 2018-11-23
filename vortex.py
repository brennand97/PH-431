import numpy as np
import matplotlib.pyplot as plt
from scipy import special as sp

# Physical constants
mu = 1      # Vacuum permeability
e0 = 1      # Vacuum permittivity
c = 1       # Speed of light

# Beam parameters
w = 1       # Beam radius
z = 1       # Z coordinate... Not used anywhere in calculations?
l = 4       # Azimuthal Laguerre index
p = 2       # Radial Laguerre index
A0 = 1      # ??

def rho(x,y):
    """
    Return the length of a vector defined on a 2D cartesian plane (like xi + yj).

    ### Parameters:
    x, y: float
    Vector components.

    ### Return Type:
    r: float
    Length of vector.
    """
    return np.sqrt(x**2+y**2) 

def u(x, y, l, p):
    """
    Evaluates the amplitude of a helical beam, at coordinates x, y relative to the beam's axis using equation 1 from Sundbeck; a special case of "Laguerre-Gaussian eigenmodes of the Helmholtz equation."

    Parameters:
    x, y: float
    Cartesian coordinates in a 2D plane, where the optical axis normal to the plane.

    l: int
    Azimuthal index of the generalized Laguerre polynomial. Also called the helicity, or topological charge of the beam.

    p: int
    Radial index of the generalized Laguerre polynomial.

    Returns:
    u: float
    Real-valued amplitude of a helical beam at coordinates x, y.

    References:
    Steven Sundbeck, Ilya Gruzberg, and David G. Grier, "Structure and scaling of helical modes of light," Opt. Lett. 30, 477-479 (2005)
    """

    # Helical beam has a radially symmetrical amplitude,
    # so the amplitude function is only dependent on the
    # distance from the origin to the x, y coordinates.
    r = rho(x,y)

    # Evaluate the equation from Sundbeck.
    return (-1)**p * (np.sqrt(2) * r/w)**l * \
        sp.genlaguerre(p, l)(2 * r**2 / w**2) * \
        np.exp(- r**2 / w**2)

def I(x, y, l, p): 
    """
    Evaluates the electromagnetic intensity of a helical beam at coordinates x, y relative to the beam's axis.

    Parameters:
    x, y: float
    Cartesian coordinates in a 2D plane, where the optical axis normal to the plane.

    l: int
    Azimuthal index of the generalized Laguerre polynomial. Also called the helicity, or topological charge of the beam.

    p: int
    Radial index of the generalized Laguerre polynomial.

    Returns:
    I: float
    Electromagnetic intensity of a helical beam at coordinates x, y.
    """
    
    return 0.5 / (mu * c) * A0**2 * ( u(x, y, l, p) )**2

# Construct a 2D space normal to the beam axis,
# over which we evaluate the beam intensity
x = np.linspace(-4,4,500)
y = x
X, Y = np.meshgrid(x,y)

plt.figure()
plt.pcolor(X,Y, I(X,Y,l,p), cmap='gray')

plt.title(r"Optical vortex intensity cross section z={}".format(z))

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

plt.show()