import numpy as np
import sys

class BogeyBailly:

    def __init__(self, grid, sigma):
        
        self.filter_width = sigma
        
        
        self.my_grid = grid
        self.start_idx = -5 #half-stencil width

        if (grid.is_2d()):
            self.ny = grid.ny
            self.twod = True

        #this is some kind of weird central difference
        self.cds_coeffs=[-0.001446093078167, \
                0.012396449873964, \
                -0.049303775636020, \
                0.120198310245186, \
                -0.199250131285813, \
                0.234810479761700,\
                -0.199250131285813, \
                0.120198310245186, \
                -0.049303775636020, \
                0.012396449873964, \
                -0.001446093078167] 



    def filter(self, f):
        
        #always filter with the cds-filter
        new_f = self.central_filter(f)

        """
        #for Euler eq. also filter with the shock sensor
        if (self.eq_type.type == "Euler"):
            new_u = self.shock_sensor(new_u, grid)
        """
        return new_f

    def shock_sensor(self, f):
        """
        not done
        """
        return f    


    def central_filter(self, f):

        #new_f = np.zeros(f.shape)
    
        #x-direction
         

        if (self.my_grid.is_1d()):
            """
            laplacian = self.spatial_operator.d_dxdx(u)
            multiplier = self.my_grid.dx * self.my_grid.dx / (np.pi*np.pi)
            u_new = u + self.filter_width * multiplier * laplacian
            """
            nx = self.my_grid.nx
            ngc = self.my_grid.ngc
            filtered = self.generic_difference(f, nx, ngc)

            return f - self.filter_width * filtered



        if(self.my_grid.is_2d()):

            filtered = np.zeros(f.shape)

            nx = self.my_grid.nx
            ny = self.my_grid.ny
            ngc = self.my_grid.ngc

            #filter x
            for j in range(ngc, ngc + ny):
                
                filtered[:, j] = self.generic_difference(f[:,j], nx, ngc )


            #filter y
            for i in range(ngc, ngc + nx):
                filtered[i, :] = self.generic_difference(f[i,:], ny, ngc )


            return f - self.filter_width * filtered

            
        
        



    def generic_difference(self, f, n, ngc):

        ret = np.zeros(f.shape)

        for i in range(len(self.cds_coeffs)):
                
                start = ngc + self.start_idx + i 
                end = start + n
                ret[ngc:-ngc] += f[start:end] * self.cds_coeffs[i]
        
        return ret
        
        
                

        



        

                



