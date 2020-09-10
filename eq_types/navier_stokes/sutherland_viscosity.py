import numpy as np



class SutherlandViscosity():

    def __init__(self, mu0, T0, C, R):

        #Sutherland coeffs (for air mu0 = 1.716E-5, T0 = 273.1, C=110.5)
        self.mu0 = mu0 
        self.T0 = T0
        self.C = C

        self.R = R

    def mu(self, T):
        
        ret = self.mu0 * (self.T0 + self.C)/(T + self.C) * (T/self.T0)**(3./2.)
        return ret

    def dmu_dx(self, diff, mu):
        
        return diff.d_dx(mu)


    def dmu_dy(self, diff, mu):

        return diff.d_dy(mu)    



    def temperature(self, p, rho):
        
        T = (p/rho)/self.R + self.T0
        return T

       
