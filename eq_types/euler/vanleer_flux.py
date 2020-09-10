import numpy as np
import sys



class VanLeer:

    
    def __init__(self, gamma):
        self.gamma = gamma
    
    
    def flux(self, U, eig, direction):
    
        fl = np.zeros(U.shape)
        fr = np.zeros(U.shape)
    
        #1d
        if (len(U.shape) == 2):
            """
            1-d conservative variables
            W = [rho, p, u]
            U = [rho, rho*E, rho*u]
            """
           
            rho = U[:,0]
            E = U[:,1] / rho
            u = U[:,2] / rho
            p = (self.gamma - 1) * rho * ( E - 0.5 * u*u)	

            a = np.sqrt(self.gamma*p/rho)
            M = u/a
            
            
            
            fl[:,0] = 0.25 * rho * a * (M + 1.0) * (M + 1.0)
            fl[:,1] = fl[:,0] * 2.0 * a * a * (1.0 + 0.5*(self.gamma-1.0)*M)**2/(self.gamma*self.gamma-1.0)
            fl[:,2] = fl[:,0] * 2.0 * a * (1.0 + 0.5*(self.gamma-1.0)*M)/self.gamma                          
            
            fr[:,0] = -0.25 * rho * a * (M - 1.0) * (M - 1.0)
            fr[:,1] = fr[:,0] * 2.0 * a * a * (1.0 - 0.5*(self.gamma-1.0)*M)**2/(self.gamma*self.gamma-1.0)
            fr[:,2] = fr[:,0] * 2.0 * a * (-1.0 + 0.5*(self.gamma-1.0)*M)/self.gamma 
    
    
            #print fr[:,2]
            #sys.exit()
    
        #2d
        if (len(U.shape) == 3):
        
            """
            2-d conservative variables
            W = [rho, p, u, v]
            U = [rho, rho*E, rho*u, rho*v]
            """
            
           

            
            rho = U[:,:,0]
            E = U[:,:,1] / rho
            u = U[:,:,2] / rho
            v = U[:,:,3] / rho
            p = (self.gamma - 1) * rho * ( E - 0.5 * (u*u + v*v))	
            

            a = np.sqrt(self.gamma*p/rho)
            
            if (direction == "x"):
                nx = 1.0; ny = 0.0
                
            if (direction == "y"):
                nx = 0.0; ny = 1.0    
            
            Vn = u*nx + v*ny
            q = np.sqrt(u*u + v*v)
            M = Vn/a
            
            temp_plus   =  (rho * a * (M + 1.0)**2)/4.0 
            temp_minus  = -(rho * a * (M - 1.0)**2)/4.0 
            
            #taken from:
            #http://www.engmech.cz/improc/2009/Bublik-201-PT.pdf
            
            gamma = self.gamma
            
            fl[:,:,0] = 1.0*temp_plus 
            fl[:,:,1] = (0.5 * (q*q - Vn) + ((gamma - 1.0)*Vn + 2*a)**2 / (2*(gamma**2 - 1.0))) * temp_plus
            fl[:,:,2] = (u + nx*(-Vn + 2*a)/gamma) * temp_plus                         
            fl[:,:,3] = (v + ny*(-Vn + 2*a)/gamma) * temp_plus 
            
            
            
            
            fr[:,:,0] = 1.0*temp_minus
            fr[:,:,1] = (0.5 * (q*q - Vn) + ((gamma - 1.0)*Vn - 2*a)**2 / (2*(gamma**2 - 1.0))) * temp_minus
            fr[:,:,2] = (u + nx*(-Vn - 2*a)/gamma) * temp_minus 
            fr[:,:,3] = (v + ny*(-Vn - 2*a)/gamma) * temp_minus
    
            #print fr[:,:,2]
            #sys.exit()
        
        return fl, fr
        
        
   
