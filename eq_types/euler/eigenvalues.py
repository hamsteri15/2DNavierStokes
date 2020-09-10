import numpy as np
 
    
def eigenvalues(gamma, W, direction):
    
    

    if (len(W.shape) == 2):
        """
        1d eigenvalues
        W = [rho, p, u]
        eig = [u-c, u, u+c]
        """
        eig = np.zeros(W.shape)
        rho = W[:,0]
        p = W[:,1]
        u = W[:,2]
        c = np.sqrt(gamma * p/rho)
        eig[:,0] = u-c
        eig[:,1] = u
        eig[:,2] = u+c
        return eig
        
    if (len(W.shape) == 3):    
        
        """
        2d eigenvalues
        W = [rho, p, u, v]
        eig_x = [u-c, u, u, u+c]
        eig_y = [v-c, v, v, v+c]
        """
        rho = W[:,:,0]
        p = W[:,:,1]
        u = W[:,:,2]
        v = W[:,:,3]
        c = np.sqrt(gamma * p/rho)
        eig = np.zeros(W.shape)

        if (direction == "x"):
            
            eig[:,:,0] = u-c
            eig[:,:,1] = u
            eig[:,:,2] = u
            eig[:,:,3] = u+c
        
        if (direction == "y"):
        
            eig[:,:,0] = v-c
            eig[:,:,1] = v
            eig[:,:,2] = v
            eig[:,:,3] = v+c
        
        
        return eig
            
            
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
