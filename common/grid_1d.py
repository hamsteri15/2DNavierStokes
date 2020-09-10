import numpy as np
import pylab as pl
import twoway_grid_stretching as stretch

class Grid:

    def __init__(self, Lx, nx, ngc, Cx=0.):
        self.nx = nx
        self.NX = nx + 2*ngc
        self.ngc = ngc
        self.num_spatial_dims = 1
        
        #domain dimensions
        self.Lx = Lx
    

        #grid is stretched
        if (Cx > 0.):
            self.stretched_x = True
            self.dx = 1./nx
        else:
            #not strect
            self.stretched_x = False
            self.dx = Lx/nx

        #stretching parameters
        self.Cx = Cx
        

        #grade the cells towards this point, set to the center for now
        self.center_x = 0.5*Lx
        

        #these are computed below in build_points()
        self.X = None
        

        #Grid first and second derivatives, computed in build_grid_derivatives()
        self.xi_x = None
        self.xi_xx = None


        
        self.build_points()
        self.build_grid_derivatives()


    def build_points(self):

        if (self.stretched_x):
            #x-dir
            xi = np.linspace(0, 1.0, self.nx)
            x = stretch.coords(xi, self.Cx, self.center_x,  self.Lx)
            
        else:
            x = np.linspace(0, self.Lx, self.nx)
        
        #physical coordinates
        self.X = x



    def build_grid_derivatives(self):

        if (self.stretched_x):
            #x-dir
            self.xi_x = stretch.first_derivative(self.X, self.Cx, self.center_x, self.Lx)
            self.xi_xx = stretch.second_derivative(self.X, self.Cx, self.center_x, self.Lx)

        else:
            self.xi_x = np.ones(self.nx)
            self.xi_xx = np.zeros(self.nx)

        #grid jacobian
        self.J = self.xi_x * 1.0



    def mult_by_xi_x(self, array):

        new_array = np.zeros(array.shape)
        ngc = self.ngc
        if (len(array.shape) == 2):
            """
            vector field
            """
            nq = array.shape[1]
            for q in range(nq):
                new_array[ngc:-ngc,q] = array[ngc:-ngc,q] * self.xi_x
                
        else:
            """
            scalar field
            """
            new_array[ngc:-ngc] = array[ngc:-ngc] * self.xi_x            

        return new_array        



    def set_grid_derivative(self, derivative):
        self.Deta_Dx = derivative 
        self.streched = True
        
    def set_point_coords(self, coords):
        self.coords = coords
    
    def is_streched(self):
        return self.streched
        
    def is_1d(self):
        return True    
    
    def is_2d(self):
        return False
    
    
    
    
     
