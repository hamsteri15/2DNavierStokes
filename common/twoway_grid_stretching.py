import numpy as np

"""
Mapping functions for the two-way grid streching. The mapping is always from 0 to 1. The functions here are taken from Anderson. EQ. X
params:
xi = xi(x) = points on an evenly distributed grid from 0 to 1
x = x(xi)    unevenly distributed points, use the function 'coords' to compute then
tau =        parameter that defines the magintude of the stretching
center =     a point (between 0 and L) in the real domain coordinates towards which the dx gets smaller
L =          real domain length  
"""


def sinh_sq(x):
    """
    sinh^2(x)
    """
    return 0.5*(np.cosh(2*x) - 1)

def coeff_B(tau, center, L):
    """
    coefficient in the formulas below
    """
    numerator = 1.0 + (np.exp(tau) - 1.0) * (center/L) 
    denominator = 1.0 + (np.exp(-tau) - 1.0) * (center/L)
    B = (1./(2*tau)) * np.log(numerator/denominator) 
    return B

#uneven coordinates (x(xi) for example)
def coords(xi, tau, center,  L):
    B = coeff_B(tau, center, L)
    return center * (1.0 + np.sinh(tau*(xi-B))/np.sinh(tau*B))

#mapping function (xi(x) for example). This is not necessary (only the derivatives are) but nice to keep here
def mapping(x, tau, center, L):
    B = coeff_B(tau, center, L)
    return B + (1./tau) * np.arcsinh( (x/center - 1.0) * np.sinh(tau*B)) 

#first derivative of the grid mapping (dxi/dx for example)
def first_derivative(x, tau, center, L):
    B = coeff_B(tau, center, L)
    numerator = np.sinh(tau * B)
    denominator = tau * center * np.sqrt(1.0 + ((x/center) - 1.0)**2 * sinh_sq(tau*B)) 

    return numerator/denominator

#second derivative of the grid mapping (ddxi/dxdx for example)
def second_derivative(x, tau, center, L):
    B = coeff_B(tau, center, L)
    return -0.1e1 / tau / center ** 2 * np.sinh(tau * B) ** 3 * ((x / center - 1) ** 2 * np.sinh(tau * B) ** 2 + 1) ** (-0.3e1 / 0.2e1) * (x / center - 1)



