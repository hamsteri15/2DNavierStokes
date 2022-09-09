import numpy as np
import sys #exit

class CDN:

    def __init__(self, grid, order):
        self.grid = grid
        self.order = order
        self.create_coeffs()

        self.name = "CDN-{0}".format(self.order)

    def create_coeffs(self):

        if (self.order == 2):
            self.d_coeffs = [-0.5, 0., 0.5]     
            self.dd_coeffs = [1.,-2.,1.]
            self.start_idx = -1
            

        elif (self.order == 4):
            self.d_coeffs = [1./12., 	-2./3., 	0., 	2./3., 	-1./12.]     
            self.dd_coeffs = [-1./12., 	4./3., 	-5./2., 	4./3., 	-1./12.]
            self.start_idx = -2
            

        elif (self.order == 6):
            self.d_coeffs = [-1./60., 	3./20., 	-3./4., 	0., 	    3./4., 	-3./20., 	1./60. ]     
            self.dd_coeffs = [1./90., 	-3./20., 	 3./2.,     -49./18., 	3./2., 	-3./20.,	1./90.]
            self.start_idx = -3
            
            
            
        else:
            print("Invalid order in CDN. Only orders 2, 4 or 6 available")




    def generic_difference(self, f, coeffs, n, ngc):

        ret = np.zeros(f.shape)

        for i in range(len(coeffs)):
                
                start = ngc + self.start_idx + i 
                end = start + n
                ret[ngc:-ngc] += f[start:end] * coeffs[i] 
        return ret
        
        
    def d_dx(self, f):
        """
        f is the function to differentiate 
        """
        if (self.grid.is_1d()):
            
            nx = self.grid.nx
            ngc = self.grid.ngc
            derivative = self.generic_difference(f, self.d_coeffs, nx, ngc)
            derivative/=self.grid.dx
                  
        if (self.grid.is_2d()):

            nx = self.grid.nx 
            ny = self.grid.ny
            ngc = self.grid.ngc
            derivative = np.zeros(f.shape)
            
            for j in range(ngc, ngc + ny):
                derivative[:, j] = self.generic_difference(f[:,j], self.d_coeffs, nx, ngc )

            derivative/=self.grid.dx

        if (self.grid.stretched_x):
            derivative = self.grid.mult_by_xi_x(derivative)
            
            
        return derivative

            
    def d_dxdx(self, f):

        if (self.grid.is_1d()):
            
            nx = self.grid.nx
            ngc = self.grid.ngc            
            #not implemented yet
            if (self.grid.stretched_x):
                df_dx = self.generic_difference(f, self.d_coeffs, nx, ngc)/self.grid.dx
                df_dxdx = self.generic_difference(f, self.dd_coeffs, nx, ngc)/(self.grid.dx*self.grid.dx)
                derivative = self.grid.mult_by_xi_xx(df_dx) + self.grid.mult_by_xi_x_sq(df_dxdx)
            else:
                derivative = self.generic_difference(f, self.dd_coeffs, nx, ngc)
                derivative/=(self.grid.dx*self.grid.dx)
            
            return derivative
                
        if (self.grid.is_2d()):
            
            nx = self.grid.nx 
            ny = self.grid.ny
            ngc = self.grid.ngc
            
            if (self.grid.stretched_x):
                df_dx = np.zeros(f.shape)
                df_dxdx = np.zeros(f.shape)
                for j in range(ngc, ngc + ny): 
                    df_dx[:, j] = self.generic_difference(f[:,j], self.d_coeffs, nx, ngc )
                    df_dxdx[:, j] = self.generic_difference(f[:,j], self.dd_coeffs, nx, ngc )

                df_dx/=self.grid.dx
                df_dxdx/=(self.grid.dx*self.grid.dx)
                derivative = self.grid.mult_by_xi_xx(df_dx) + self.grid.mult_by_xi_x_sq(df_dxdx)
                

            else:
                derivative = np.zeros(f.shape)
                for j in range(ngc, ngc + ny): 
                    derivative[:, j] = self.generic_difference(f[:,j], self.dd_coeffs, nx, ngc )

                derivative/=(self.grid.dx*self.grid.dx)

            return derivative
            
                      
    
    
    def d_dy(self, f):
        
        if (self.grid.is_1d()):
            print("Trying to compute d_dy on a 1D grid. Dont!")

            sys.exit()

        nx = self.grid.nx 
        ny = self.grid.ny
        ngc = self.grid.ngc
        derivative = np.zeros(f.shape)
        
        for i in range(ngc, ngc + nx):
            derivative[i,:] = self.generic_difference(f[i,:], self.d_coeffs, ny, ngc )

        derivative/=(self.grid.dy)

        
        if (self.grid.stretched_y):
            derivative = self.grid.mult_by_eta_y(derivative)
        
        return derivative

        



    def d_dydy(self, f):
        
        if (self.grid.is_1d()):
            print ("Trying to compute d_dydy on a 1D grid. Dont!")
            sys.exit()

        nx = self.grid.nx 
        ny = self.grid.ny
        ngc = self.grid.ngc

        if (self.grid.stretched_y):
            df_dy = np.zeros(f.shape)
            df_dydy = np.zeros(f.shape)
            for i in range(ngc, ngc + nx):
                df_dy[i,:] = self.generic_difference(f[i,:], self.d_coeffs, ny, ngc )
                df_dydy[i,:] = self.generic_difference(f[i,:], self.dd_coeffs, ny, ngc )

            df_dy/=self.grid.dy
            df_dydy/=(self.grid.dy*self.grid.dy)
            derivative = self.grid.mult_by_eta_yy(df_dy) + self.grid.mult_by_eta_y_sq(df_dydy)
            

        else:

            derivative = np.zeros(f.shape)
            for i in range(ngc, ngc + nx):
                derivative[i,:] = self.generic_difference(f[i,:], self.dd_coeffs, ny, ngc )

            derivative/=(self.grid.dy*self.grid.dy)

            
        return derivative
        
        
        


                        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
                
