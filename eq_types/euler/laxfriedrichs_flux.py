import numpy as np
from physical_flux import PhysicalFlux



class LaxFriedrichs:

    
    def __init__(self, gamma):
        self.gamma = gamma
    
    
    def flux(self, U, eig, direction):
        """
        Computes the left and right split fluxes
        """
    
        physical_flux = PhysicalFlux(self.gamma)

        F = physical_flux.flux(U, direction)        
        
       
        #this is the maximum of the whole domain
        alpha = np.amax(np.abs(eig))

        fl = 0.5*(F + alpha*U)
        fr = 0.5*(F - alpha*U)
            
        return fl, fr 

    
        
    

    

    
    
        
        
