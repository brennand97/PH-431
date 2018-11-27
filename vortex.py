import numpy as np
import matplotlib.pyplot as plt
from scipy import special as sp

# Physical constants
mu = 1      # Vacuum permeability
e0 = 1      # Vacuum permittivity
c = 1       # Speed of light

# Beam parameters
w = 1       # Beam radius
z = 1       # Z coordinate
t = 0       # time
l = 4       # Azimuthal Laguerre index
p = 2       # Radial Laguerre index
A0 = 1      # Initial magnitude
k = 1       # wave vector
omega = 1   # frequency

E_th = 0
E_u = [np.cos(E_th), np.sin(E_th)]
B_u = [np.cos(E_th + np.pi / 2), np.sin(E_th + np.pi / 2)]

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
    r = rho(x,y)

    # Evaluate the equation from Sundbeck.
    return (-1)**p * (np.sqrt(2) * r/w)**l * \
        sp.genlaguerre(p, l)(2 * r**2 / w**2) * \
        np.exp(- r**2 / w**2)

def I(x, y, l, p): 
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
    
    return 0.5 / (mu * c) * A0**2 * ( u(x, y, l, p) )**2

def S(x, y, l, p):
    return e0 * c * A0**2 * ( u(x, y, l, p) )**2 * np.cos(np.arctan2(x, y) * l + k*z - omega*t)**2

def E_mag(x, y, l, p):
    return A0 * u(x, y, l, p) * np.cos(np.arctan2(x, y) * l + k*z - omega*t)

def B_mag(x, y, l, p):
    return A0 / c * u(x, y, l, p) * np.cos(np.arctan2(x, y) * l + k*z - omega*t)

def E(x, y, l, p):
    mag = E_mag(x, y, l, p)
    return E_u[0] * mag, E_u[1] * mag

def B(x, y, l, p):
    mag = B_mag(x, y, l, p)
    return B_u[0] * mag, B_u[1] * mag

# Construct a 2D space normal to the beam axis,
# over which we evaluate the beam intensity
x = np.linspace(-4,4,1000)
y = x
X, Y = np.meshgrid(x,y)



plt.figure()
plt.pcolor(X,Y, I(X,Y,l,p), cmap='hot')

plt.colorbar(label = "Intensity")
plt.title("Optical vortex intensity cross section z={}".format(z))
plt.xlabel('x')
plt.ylabel('y')

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

plt.savefig('intensity.png', dpi=500)



plt.figure()
plt.pcolor(X,Y, S(X,Y,l,p), cmap='hot')

plt.colorbar(label = "Magnitude")
plt.title("Optical vortex poynting vector z={}".format(z))
plt.xlabel('x')
plt.ylabel('y')

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

plt.savefig('poynting.png', dpi=500)



plt.figure()
plt.pcolor(X,Y, E_mag(X,Y,l,p), cmap='RdBu')

plt.colorbar(label = "Magnitude")
plt.title("Optical vortex electric magnitude z={}".format(z))
plt.xlabel('x')
plt.ylabel('y')

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

plt.savefig('e_mag.png', dpi=500)



plt.figure()
plt.pcolor(X,Y, B_mag(X,Y,l,p), cmap='RdBu')

plt.colorbar(label = "Magnitude")
plt.title("Optical vortex magnetic magnitude z={}".format(z))
plt.xlabel('x')
plt.ylabel('y')

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

plt.savefig('b_mag.png', dpi=500)



x = np.linspace(-4,4,200)
y = x
X, Y = np.meshgrid(x,y)



plt.figure()
plt.quiver(X,Y, *E(X,Y,l,p))

plt.title("Optical vortex electric field z={}".format(z))
plt.xlabel('x')
plt.ylabel('y')

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

plt.savefig('e_field.png', dpi=500)



plt.figure()
plt.quiver(X,Y, *B(X,Y,l,p))

plt.title("Optical vortex magnetic field z={}".format(z))
plt.xlabel('x')
plt.ylabel('y')

ax = plt.gca()
ax.set_aspect('equal', adjustable='box')

plt.savefig('b_field.png', dpi=500)

plt.show()