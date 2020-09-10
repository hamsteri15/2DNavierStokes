import numpy as np
import pylab as pl
import twoway_grid_stretching as stretch


class Grid:

    

    def __init__(self, params):

        #cell counts with and without gcs and and al
        self.nx = params.nx
        self.NX = params.nx + 2*params.ngc
        self.ny = params.ny
        self.NY = params.ny + 2*params.ngc
        self.ngc = params.ngc #number of ghost points

        self.num_spatial_dims = 2
        
        #domain dimensions
        self.Lx = params.Lx
        self.Ly = params.Ly
        
        #stretching parameters
        self.Cx = params.Cx
        self.Cy = params.Cy 



        #stepsizes on the even grid mapping always from 0 to 1 so the length 
        #in the computational space is 1

        if (self.Cx > 0.):
            self.stretched_x = True
            self.dx = 1./self.nx
        else:
            self.stretched_x = False
            self.dx = self.Lx/self.nx


        if (self.Cy > 0.):
            self.stretched_y = True
            self.dy = 1./self.ny
        else:
            self.stretched_y = False
            self.dy = self.Ly/self.ny


        
        
        

        #grade the cells towards this point, set to the center for now
        self.center_x = 0.5*self.Lx
        self.center_y = 0.5*self.Ly

        #these are computed below in build_points()
        self.X = None
        self.Y = None

        #Grid first and second derivatives, computed in build_grid_derivatives()
        self.xi_x = None
        self.eta_y = None
        self.xi_xx = None
        self.eta_yy = None

        
        self.build_points()
        self.build_grid_derivatives()
        

    def build_points(self):

        if (self.stretched_x):
            #x-dir
            xi = np.linspace(0, 1.0, self.nx)
            x = stretch.coords(xi, self.Cx, self.center_x,  self.Lx)
            
        else:
            x = np.linspace(0, self.Lx, self.nx)

        if (self.stretched_y):
            #y-dir
            eta = np.linspace(0, 1.0, self.ny)
            y = stretch.coords(eta, self.Cy, self.center_y,  self.Ly)
            
        else:
            y = np.linspace(0, self.Ly, self.ny)

        #physical coordinates
        self.X,self.Y = np.meshgrid(x,y)
        
       
        


    def build_grid_derivatives(self):

        if (self.stretched_x):
            #x-dir
            self.xi_x = stretch.first_derivative(self.X, self.Cx, self.center_x, self.Lx)
            self.xi_xx = stretch.second_derivative(self.X, self.Cx, self.center_x, self.Lx)

        else:
            self.xi_x = np.ones((self.nx, self.ny))
            self.xi_xx = np.zeros((self.nx, self.ny))


        if (self.stretched_y):
            #y-dir
            self.eta_y = stretch.first_derivative(self.Y, self.Cy, self.center_y, self.Ly)
            self.eta_yy = stretch.second_derivative(self.Y, self.Cy, self.center_y, self.Ly)

        else:
            self.eta_y = np.ones((self.nx, self.ny))
            self.eta_yy = np.zeros((self.nx, self.ny))

        
        #transpose to get the correct order of rows and columns
        self.xi_x = np.transpose(self.xi_x)
        self.xi_xx = np.transpose(self.xi_xx)
        self.eta_y = np.transpose(self.eta_y)
        self.eta_yy = np.transpose(self.eta_yy)
        self.X = np.transpose(self.X)
        self.Y = np.transpose(self.Y)
        self.J = self.xi_x * self.eta_y

    def mult_by_xi_x(self, array):

        new_array = np.zeros(array.shape)
        ngc = self.ngc
        
        if (len(array.shape) == 3):
            nq = array.shape[2]
            for q in range(nq):
                new_array[ngc:-ngc, ngc:-ngc, q] = array[ngc:-ngc, ngc:-ngc, q] * self.xi_x
                  
        else:
                
            new_array[ngc:-ngc, ngc:-ngc] = array[ngc:-ngc, ngc:-ngc] * self.xi_x
             
        return new_array

    def mult_by_eta_y(self, array):

        new_array = np.zeros(array.shape)
        ngc = self.ngc
        if (len(array.shape) == 3):
            nq = array.shape[2]
            for q in range(nq):    
                new_array[ngc:-ngc, ngc:-ngc, q] = array[ngc:-ngc, ngc:-ngc, q] * self.eta_y
    
        else:
                 
            new_array[ngc:-ngc, ngc:-ngc] = array[ngc:-ngc, ngc:-ngc] * self.eta_y

        return new_array  

    def mult_by_xi_xx(self, array):
        """
        multiplies given input array by xi_xx, usually called for second derivatives for which
        df_dxdx = xi_xx * df_dxi + xi_x*xi_x * df_dxixi 
        """
        new_array = np.zeros(array.shape)
        ngc = self.ngc
        if (len(array.shape) == 3):
            nq = array.shape[2]
            for q in range(nq):
                
                new_array[ngc:-ngc, ngc:-ngc, q] = array[ngc:-ngc, ngc:-ngc, q] * self.xi_xx
        else:
            new_array[ngc:-ngc, ngc:-ngc] = array[ngc:-ngc, ngc:-ngc] * self.xi_xx


        return new_array

    def mult_by_xi_x_sq(self, array):
        """
        multiplies given input array by xi_x*xi_x, usually called for second derivatives for which
        df_dxdx = xi_xx * df_dxi + xi_x*xi_x * df_dxixi 
        """
        new_array = np.zeros(array.shape)
        ngc = self.ngc
        if (len(array.shape) == 3):
            nq = array.shape[2]
            for q in range(nq):
                
                new_array[ngc:-ngc, ngc:-ngc, q] = array[ngc:-ngc, ngc:-ngc, q] * self.xi_x*self.xi_x
        else:
            new_array[ngc:-ngc, ngc:-ngc] = array[ngc:-ngc, ngc:-ngc] * self.xi_x*self.xi_x


        return new_array

    def mult_by_eta_yy(self, array):
        
        new_array = np.zeros(array.shape)
        ngc = self.ngc
        if (len(array.shape) == 3):
            nq = array.shape[2]
            for q in range(nq):
                
                new_array[ngc:-ngc, ngc:-ngc, q] = array[ngc:-ngc, ngc:-ngc, q] * self.eta_yy
        else:
            new_array[ngc:-ngc, ngc:-ngc] = array[ngc:-ngc, ngc:-ngc] * self.eta_yy

        return new_array

    def mult_by_eta_y_sq(self, array):
        
        new_array = np.zeros(array.shape)
        ngc = self.ngc
        if (len(array.shape) == 3):
            nq = array.shape[2]
            for q in range(nq):
                
                new_array[ngc:-ngc, ngc:-ngc, q] = array[ngc:-ngc, ngc:-ngc, q] * self.eta_y*self.eta_y
        else:
            new_array[ngc:-ngc, ngc:-ngc] = array[ngc:-ngc, ngc:-ngc] * self.eta_y*self.eta_y

        return new_array



    def is_1d(self):
        return False    
    
    def is_2d(self):
        return True
        
        
        








    #######################
    #Visualization routines
    #######################    
        
        

    def visualize(self):

        # Creates two subplots and unpacks the output array immediately
        f, axes = pl.subplots(3, 2, figsize=(10,15)) #rows, cols

        self.plot_grid(axes[0,0])
        self.plot_jacobian(axes[0,1])
        self.plot_xix(axes[1,0])
        self.plot_etay(axes[1,1])
        #self.plot_xix_per_J(axes[2,0])
        #self.plot_etay_per_J(axes[2,1])
        self.plot_xixx(axes[2,0])
        self.plot_etayy(axes[2,1])
        pl.show()
        
    
    


    
        
    
    
    



    def plot_grid(self, ax):

        temp = np.zeros(self.X.shape) 
        im = ax.pcolor(self.X,self.Y,temp, edgecolors="red" )
        pl.colorbar(im,ax=ax)
        ax.set_title("Gridlines")
        ax.set_xlim(np.amin(self.X), np.amax(self.X))
        ax.set_ylim(np.amin(self.Y), np.amax(self.Y))

    def plot_jacobian(self, ax):
        im = ax.pcolor(self.X,self.Y,self.J )
        pl.colorbar(im,ax=ax)
        ax.set_title("Jacobian")
        ax.set_xlim(np.amin(self.X), np.amax(self.X))
        ax.set_ylim(np.amin(self.Y), np.amax(self.Y))

    def plot_xix(self, ax):
        im = ax.pcolor(self.X,self.Y,self.xi_x)
        pl.colorbar(im,ax=ax)
        ax.set_title("$\\xi_x$")
        ax.set_xlim(np.amin(self.X), np.amax(self.X))
        ax.set_ylim(np.amin(self.Y), np.amax(self.Y))

    def plot_xixx(self, ax):
        im = ax.pcolor(self.X,self.Y,self.xi_xx)
        pl.colorbar(im,ax=ax)
        ax.set_title("$\\xi_{xx}$")
        ax.set_xlim(np.amin(self.X), np.amax(self.X))
        ax.set_ylim(np.amin(self.Y), np.amax(self.Y))

    def plot_etay(self, ax):
        im = ax.pcolor(self.X,self.Y,self.eta_y)
        pl.colorbar(im,ax=ax)
        ax.set_title("$\\eta_y$")
        ax.set_xlim(np.amin(self.X), np.amax(self.X))
        ax.set_ylim(np.amin(self.Y), np.amax(self.Y))


    def plot_etayy(self, ax):
        im = ax.pcolor(self.X,self.Y,self.eta_yy)
        pl.colorbar(im,ax=ax)
        ax.set_title("$\\eta_{yy}$")
        ax.set_xlim(np.amin(self.X), np.amax(self.X))
        ax.set_ylim(np.amin(self.Y), np.amax(self.Y))

    def plot_xix_per_J(self, ax):
        im = ax.pcolor(self.X,self.Y,self.xi_x/self.J)
        pl.colorbar(im,ax=ax)
        ax.set_title("$\\xi_x/J$")
        ax.set_xlim(np.amin(self.X), np.amax(self.X))
        ax.set_ylim(np.amin(self.Y), np.amax(self.Y))


    def plot_etay_per_J(self, ax):
        im = ax.pcolor(self.X,self.Y,self.eta_y/self.J)
        pl.colorbar(im,ax=ax)
        ax.set_title("$\\eta_y/J$")
        ax.set_xlim(np.amin(self.X), np.amax(self.X))
        ax.set_ylim(np.amin(self.Y), np.amax(self.Y))
