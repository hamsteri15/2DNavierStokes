import numpy as np
import sys

from eigenvalues import eigenvalues
from grid_2d import Grid


class Euler:

    def __init__(self, grid_parameters, time_parameters, solver_parameters, fluid_properties):


        self.type = "Euler"
        self.grid = Grid(grid_parameters)
        
        self.gamma = fluid_properties.gamma
        self.time_scheme = time_parameters.get_time_scheme()
        self.spatial_operator_convection = solver_parameters.get_convection_scheme(self.grid)
        self.flux_method = solver_parameters.get_convection_flux(self.gamma)
        self.dt = time_parameters.dt
        self.adjust_dt = time_parameters.adjust_dt
        
        self.filter_method = solver_parameters.get_filter(self.grid)
        
        #eigenvalues in x- and y-directions
        self.eig_x = None
        self.eig_y = None 

        
        self.d_dt = self.time_scheme.d_dt

        if (self.grid.is_1d()):
            self.oned = True
            self.dU = self.dU_convection_1d
        else:
            self.oned = False    
            self.dU = self.dU_convection_2d
            
   

    def dU_convection_1d(self, U):
        """
        This function specifies how the Euler equations are differentiated. Whether to differentiate the physical 
        Euler flux directly, transform into characteristic space and differentiate there and transform back, use a Riemann solver
        etc. 
        """

        d_dx = self.spatial_operator_convection.d_dx


        if ("Weno" in self.spatial_operator_convection.name):
            """
            Weno scheme always solves for the split fluxes
            """
            
            #flux_type = EulerFluxSplit(self.gamma, riemann_solver)
            fl, fr = self.flux_method.flux(U, self.eig_x, direction="x")
            
            return -d_dx(fl, fr)

            


        else:
            """
            Other schemes discretize the flux vector directly
            """
            
            F = self.flux_method.flux(U, direction="x")
            return -d_dx(F)


    def dU_convection_2d(self, U):
        """
        This function specifies how the Euler equations are differentiated. Whether to differentiate the physical 
        Euler flux directly, transform into characteristic space and differentiate there and transform back, use a Riemann solver
        etc. 
        """

        d_dx = self.spatial_operator_convection.d_dx
        d_dy = self.spatial_operator_convection.d_dy


        if ("Weno" in self.spatial_operator_convection.name):
            """
            Weno scheme solves for the split fluxes
            """

            fl, fr = self.flux_method.flux(U, self.eig_x, direction="x")
            gl, gr = self.flux_method.flux(U, self.eig_y, direction="y")
            Rx = d_dx(fl, fr)
            Ry = d_dy(gl, gr)
            R = -(Rx + Ry)
            return R
            
            



        else:
            """
            Other schemes discretize the flux vector directly
            """
            
            Fx = self.flux_method.flux(U, direction="x")
            Fy = self.flux_method.flux(U, direction="y")

            return -(d_dx(Fx) + d_dy(Fy))


    def get_eigenvalues(self, W):
        
        eig_x = eigenvalues(self.gamma, W, "x")
        eig_y = eigenvalues(self.gamma, W, "y")
        return eig_x, eig_y
    


    def take_step(self, W):       

        self.eig_x, self.eig_y = self.get_eigenvalues(W)
        
        if (self.adjust_dt):
            self.dt = self.adjust_timestep(self.eig_x, self.eig_y)
            
        

        U = self.primitive_to_conservative(W)
        increment = self.d_dt(U, self.dU, self.dt)
        Unew = U + increment
        Wnew = self.conservative_to_primitive(Unew)

        Wfiltered = self.filter_method.filter(Wnew)
        return Wfiltered



       
    def adjust_timestep(self, eig_x, eig_y=None):
        """
        adjusts the timestep based on the maximum wave speeds
        """
        max_cfl = 0.4
        dx = self.grid.dx
        max_eig_x = np.amax(np.abs(eig_x))
        dt_x = dx / max_eig_x
        #print max_eig_x
        


        if (self.grid.is_2d()):

            dy = self.grid.dy            
            max_eig_y = np.amax(np.abs(eig_y))
            dt_y = dy / max_eig_y
            dt = max_cfl * min(dt_x, dt_y)
            
        else:
            dt = max_cfl * dt_x
	    

        return dt

    def primitive_to_conservative(self, W):
    
        U = np.zeros(W.shape)
        if (len(W.shape) == 2):
            """
            1-d conservative variables
            W = [rho, p, u]
            U = [rho, rho*E, rho*u]
            """
            rho = W[:,0]
            p = W[:,1]
            u = W[:,2]    
            E = p / ((self.gamma-1)*rho) + 0.5*u*u
            
            U[:,0] = rho
            U[:,1] = rho*E
            U[:,2] = rho*u
            
    
        if (len(W.shape) == 3):
            """
            2-d conservative variables
            W = [rho, p, u, v]
            U = [rho, rho*E, rho*u, rho*v]
            """
            
            U[:,:,0] = W[:,:,0]             #rho
            U[:,:,2] = W[:,:,0] * W[:,:,2]  #rho*u
            U[:,:,3] = W[:,:,0] * W[:,:,3]  #rho*v
            kin = 0.5 * W[:,:,0] * ( W[:,:,2]*W[:,:,2] + W[:,:,3]*W[:,:,3])
            U[:,:,1] = W[:,:,1] / (self.gamma - 1.0) + kin
            
            
        return U     
            
            
    
    def conservative_to_primitive(self, U):
    
    
        if (len(U.shape) == 2):
            """
            1-d conservative variables
            W = [rho, p, u]
            U = [rho, rho*E, rho*u]
            """
            W = np.zeros(U.shape)

            rho = U[:,0]
            E = U[:,1] / rho
            u = U[:,2] / rho
            p = (self.gamma - 1) * rho * ( E - 0.5 * u*u)	

            

            W[:,0] = rho
            W[:,1] = p
            W[:,2] = u

            return W
            
    
        if (len(U.shape) == 3):
            """
            2-d conservative variables
            W = [rho, p, u, v]
            U = [rho, rho*E, rho*u, rho*v]
            """
            W = np.zeros(U.shape)
            
            W[:,:,0] = U[:,:,0]             #rho
            W[:,:,2] = U[:,:,2] / U[:,:,0]  #u
            W[:,:,3] = U[:,:,3] / U[:,:,0]  #v
            
            
            kin = 0.5 * (W[:,:,2]*W[:,:,2] + W[:,:,3]*W[:,:,3] )
            W[:,:,1] = (self.gamma - 1.0) * (U[:,:,1] - U[:,:,0] * kin)            
            
            
            return W
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
