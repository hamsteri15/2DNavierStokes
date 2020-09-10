import numpy as np
import sys


class Weno5:

    

    def __init__(self, grid):
        
        self.grid = grid
        self.weight_type = "Shu"
        self.name = "Weno-5"

        
        
    

        

    def weno_left(self, fL_m2, fL_m1, fL_, fL_p1, fL_p2):


	
        #ENO-stencils
        f0 = (2.*fL_m2 - 7.*fL_m1 + 11.*fL_)*(1./6.)
        f1 = (-fL_m1 + 5.*fL_ + 2.*fL_p1)*(1./6.)
        f2 = (2.*fL_ + 5.*fL_p1 - fL_p2)*(1./6.)

        d_left  = [ 0.1, 0.6, 0.3 ]
         

        #if (weights == "Shu"):

        weights = self.weights_shu(fL_m2, fL_m1, fL_, fL_p1, fL_p2, d_left)

        #else ...

        return weights[0]*f0 + weights[1]*f1 + weights[2]*f2;

        
    def weno_right(self, fR_m2, fR_m1, fR_, fR_p1, fR_p2):


        #ENO-stencils
        f0 = (-fR_m2 + 5.*fR_m1 + 2.*fR_)*(1./6.)
        f1 = (2.*fR_m1 + 5.*fR_ - fR_p1)*(1./6.)
        f2 = (11.*fR_ - 7.*fR_p1 + 2.*fR_p2)*(1./6.)


        d_right = [ 0.3, 0.6, 0.1 ]
        weights = self.weights_shu(fR_m2, fR_m1, fR_, fR_p1, fR_p2, d_right)	


        return weights[0]*f0 + weights[1]*f1 + weights[2]*f2;
	


    def weights_shu(self, f_m2, f_m1, f_, f_p1, f_p2, d):

        
        temp1 = f_m2 - 2.*f_m1 + f_
        temp2 = f_m2 - 4.*f_m1 + 3.*f_
        b0 = (13./12.) * temp1 * temp1 + (1./4.) * temp2 * temp2

        temp1 = f_m1 - 2.*f_ + f_p1;
        temp2 = f_m1 - f_p1;
        b1 = (13./12.) * temp1 * temp1 + (1./4.) * temp2 * temp2


        temp1 = f_ - 2.*f_p1 + f_p2
        temp2 = 3.*f_ -4.*f_p1 + f_p2
        b2 = (13./12.) * temp1 * temp1 + (1./4.) * temp2 * temp2	


        epsilon_shu = 1E-6

        alpha_0 = d[0]/((epsilon_shu + b0)*(epsilon_shu + b0))
        alpha_1 = d[1]/((epsilon_shu + b1)*(epsilon_shu + b1))
        alpha_2 = d[2]/((epsilon_shu + b2)*(epsilon_shu + b2))

        alpha = alpha_0 + alpha_1 + alpha_2;


        return [alpha_0/alpha, alpha_1/alpha, alpha_2/alpha]

	

    def oned_numflux_left(self, f, grid):
        """
        Get the stencils i - 3 to i + 2 and pass to the weno routine
        """

        nx = grid.nx
        ngc = grid.ngc
        width = nx + 1 #flux stencils located at i+1/2, therefore width is nx + 1
        stencils = []
        for i in range(5):

            start_i = ngc - 3 + i
            end_i = start_i + width
            stencil = f[start_i : end_i]
            stencils.append(stencil)

        numflux = self.weno_left(stencils[0], stencils[1], stencils[2], stencils[3], stencils[4])
        return numflux
   
    def oned_numflux_right(self, f, grid):
        """
        Get the stencils i - 2 to i + 3 and pass to the weno routine
        """
        nx = grid.nx
        ngc = grid.ngc
        width = nx + 1 #flux stencils located at i+1/2, therefore width is nx + 1
        stencils = []
        for i in range(5):

            start_i = ngc - 2 + i
            end_i = start_i + width
            stencil = f[start_i : end_i]
            stencils.append(stencil)

        numflux = self.weno_right(stencils[0], stencils[1], stencils[2], stencils[3], stencils[4])
        return numflux
   
   
    def twod_numflux_left(self, f, grid, direction):
        

        if (direction == "x"):
        
            nx = grid.nx
            ngc = grid.ngc
            width = nx + 1 #flux stencils located at i+1/2, therefore width is nx + 1
            stencils = []
            for i in range(5):

                start_i = ngc - 3 + i
                end_i = start_i + width
                stencil = f[start_i : end_i, ngc:-ngc, :]
                stencils.append(stencil)


        if (direction == "y"):
        
            ny = grid.ny
            ngc = grid.ngc
            width = ny + 1 #flux stencils located at i+1/2, therefore width is nx + 1
            stencils = []
            for j in range(5):
                start_j = ngc - 3 + j
                end_j = start_j + width
                stencil = f[ngc:-ngc, start_j : end_j,  :]
                stencils.append(stencil)

        numflux = self.weno_left(stencils[0], stencils[1], stencils[2], stencils[3], stencils[4])
    
        return numflux
   
    
    def twod_numflux_right(self, f, grid, direction):
        
        if (direction == "x"):
        
            nx = grid.nx
            ngc = grid.ngc
            width = nx + 1 #flux stencils located at i+1/2, therefore width is nx + 1
            stencils = []
            for i in range(5):

                start_i = ngc - 2 + i
                end_i = start_i + width
                stencil = f[start_i : end_i, ngc:-ngc, :]
                stencils.append(stencil)
            

        if (direction == "y"):
        
            ny = grid.ny
            ngc = grid.ngc
            width = ny + 1 #flux stencils located at i+1/2, therefore width is nx + 1
            stencils = []
            for j in range(5):
                start_j = ngc - 2 + j
                end_j = start_j + width
                stencil = f[ngc:-ngc, start_j : end_j,  :]
                stencils.append(stencil)
            

            
    
    
        numflux = self.weno_right(stencils[0], stencils[1], stencils[2], stencils[3], stencils[4])
    
        return numflux
    
    

    def d_dx(self, fl, fr):
        """
        generic weno x-derivative
        """
        ngc = self.grid.ngc
        derivative = np.zeros(fl.shape)
        

        if (self.grid.is_1d()):
            
            #Left and right biased weno fluxes
            f_left = self.oned_numflux_left(fl,self.grid)
            f_right = self.oned_numflux_right(fr, self.grid)
    
            temp = f_left + f_right
        
            derivative[ngc:-ngc,:] = (temp[1:,:] - temp[0:-1,:])/self.grid.dx  
              
            
                
        if (self.grid.is_2d()):
            
            #Weno flux
            f_left = self.twod_numflux_left(fl, self.grid, "x")
            f_right = self.twod_numflux_right(fr, self.grid, "x")
            
            temp = f_left + f_right
            
            derivative[ngc:-ngc,ngc:-ngc,:] = np.diff(temp, axis=0)/self.grid.dx

        
        if (self.grid.stretched_x):
            derivative = self.grid.mult_by_xi_x(derivative)
            
        

        return derivative

    
     
            

    
    def d_dy(self, gl, gr):
        """
        generic weno y-derivative
        """

        if (self.grid.is_1d()):
            #arbitrary, never differentiate in y-direction for 1D cases
            derivative = np.zeros(u.dimensions)
            return derivative
                
                
        if (self.grid.is_2d()):
            ngc = self.grid.ngc
            derivative = np.zeros(gl.shape)

            
        
            #Weno flux
            g_left = self.twod_numflux_left(gl, self.grid, "y")
            g_right = self.twod_numflux_right(gr, self.grid, "y")
            
            temp = g_left + g_right
            
            derivative[ngc:-ngc,ngc:-ngc,:] = np.diff(temp, axis=1)/self.grid.dy

            
        if (self.grid.stretched_y):
            derivative = self.grid.mult_by_eta_y(derivative)
            
            
        return derivative     
                        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
                
