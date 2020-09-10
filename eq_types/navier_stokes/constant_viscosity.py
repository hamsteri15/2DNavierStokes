import numpy as np



class ConstantViscosity():

    def __init__(self, mu0, R):

        self.mu0 = mu0
        self.R = R #gas constant

    def mu(self, T):
        return self.mu0

    def dmu_dx(self, diff, mu):
        return 0.

    def dmu_dy(self, diff, mu):
        return 0.    

    def temperature(self, p, rho):
        
        T = (p/rho)/self.R
        return T

