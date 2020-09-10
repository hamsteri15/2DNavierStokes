import numpy as np



class PhysicalFlux:
    """
    Given conservative variables U, computes the Euler equation physical flux
    """


    def __init__(self, gamma):
        self.gamma = gamma
        

    
    def flux(self, U, direction):
        """
        Computes the physical flux of the Euler equations, some schemes differentiate this 
        directly
        """
        
        if (len(U.shape) == 2):
            """
            1-d Euler flux
            U = [rho, rho*E, rho*u]
            F(U) = [rho*u, u*(rho*E + p), rho*u*u + p]
            """
            F = np.zeros(U.shape)
            rho = U[:,0]
            rhoE = U[:,1]
            rhou = U[:,2]
            
            E = rhoE / rho
            u = rhou / rho
            p = (self.gamma - 1) * rho * ( E - 0.5 * u*u)
            
            
            F[:,0] = rhou
            F[:,1] = u*(rhoE + p)  
            F[:,2] = rhou*u + p
            return F
            
        if (len(U.shape) == 3):
            """
            2-d Euler flux
            """

            flux = np.zeros(U.shape)
            rho = U[:,:,0]
            rhoE = U[:,:,1]
            rhou = U[:,:,2]
            rhov = U[:,:,3]
            u = rhou/rho
            v = rhov/rho
            E = rhoE/rho
            p = (self.gamma - 1) * rho * ( E - 0.5 * (u*u + v*v))


            if (direction == "x"):
                flux[:,:,0] = rhou
                flux[:,:,1] = u*(rhoE +p)  
                flux[:,:,2] = rhou*u + p
                flux[:,:,3] = rhou*v

            if (direction == "y"):
                flux[:,:,0] = rhov
                flux[:,:,1] = v*(rhoE +p)  
                flux[:,:,2] = rhov*u
                flux[:,:,3] = rhov*v + p

            return flux
            
            
            
            
            
            
        

    
        
            
            
    
    
    
    
    
    
            
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
