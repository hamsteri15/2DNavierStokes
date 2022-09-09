import numpy as np
import sys

sys.path.append('../boundary_conditions')
from boundary_neumann import BoundaryNeumann
from boundary_periodic import BoundaryPeriodic


class TestCase:


    def __init__(self, name, boundary_conditions):
        
        self.name = name
        self.boundary_conditions = boundary_conditions
        
        
    def apply_boundary_conditions(self, W, grid):
    
        for boundary in self.boundary_conditions:
            boundary.apply(W, grid)


    def apply_neumann_all(self, W, grid):

        #apply neumann bcs to initial conditions to fill up the corners for 
        #cross-derivatives
        temp_boundaries = [BoundaryNeumann("left"), BoundaryNeumann("right"), BoundaryNeumann("top"), BoundaryNeumann("bottom")]
        for boundary in temp_boundaries:
            boundary.apply(W,grid)


class Sod_x(TestCase):
    def __init__(self, name):
        boundary_conditions = [BoundaryNeumann("left"), BoundaryNeumann("right"), BoundaryNeumann("top"), BoundaryNeumann("bottom")]
    
        TestCase.__init__(self, name, boundary_conditions)

        
    def get_initial_condition(self, grid):
        
        X = grid.X
        Y = grid.Y
        nx = grid.nx
        ny = grid.ny
        ngc = grid.ngc
        Wtemp = np.zeros((nx, ny, 4))

        #rho p u v
        left=  [1.0, 1.0, 0.0, 0.0]  
        right = [0.125, 0.1, 0.0, 0.0]
        Wtemp[0:int(nx/2), :, :] = left
        Wtemp[int(nx/2):, :, :] = right
        


        W = np.zeros((nx + 2*ngc, ny + 2*ngc, 4))
        W[ngc:-ngc, ngc:-ngc, :] = Wtemp
        
        self.apply_neumann_all(W, grid)
        
        return W         


class QuadrantTestCase(TestCase):
    
    def __init__(self, name, upper_left, upper_right, lower_left, lower_right):
    
        boundary_conditions = [BoundaryNeumann("left"), BoundaryNeumann("right"), BoundaryNeumann("top"), BoundaryNeumann("bottom")]
    
        TestCase.__init__(self, name, boundary_conditions)
    
        self.upper_left = upper_left
        self.upper_right = upper_right
        self.lower_left = lower_left
        self.lower_right = lower_right
    
        
    def get_initial_condition(self, grid):
        
        X = grid.X
        Y = grid.Y
        nx = grid.nx
        ny = grid.ny
        ngc = grid.ngc
        Wtemp = np.zeros((nx, ny, 4))

        Wtemp[0:int(nx/2), 0:int(ny/2), :] = self.upper_left
        Wtemp[int(nx/2):, 0:int(ny/2), :] = self.upper_right
        Wtemp[0:int(nx/2):, int(ny/2):, :] = self.lower_left    
        Wtemp[int(nx/2):, int(ny/2):, :] =  self.lower_right
        


        W = np.zeros((nx + 2*ngc, ny + 2*ngc, 4))
        W[ngc:-ngc, ngc:-ngc, :] = Wtemp
        
        self.apply_neumann_all(W, grid)
        
        return W
       
class Case3(QuadrantTestCase):

    def __init__(self):
    
        #rho p u v
        upper_left=  [0.5323, 0.3, 1.206, 0.0]  
        upper_right= [1.5, 1.5, 0.0, 0.0]
        lower_left = [0.138, 0.029, 1.206, -1.206]
        lower_right= [0.5323, 0.3, 0.0, -1.206]       
       
        QuadrantTestCase.__init__(self, "Case3", upper_left, upper_right, lower_left, lower_right)
        

class Case4(QuadrantTestCase):

    def __init__(self):
    
        #rho p u v
        upper_left=  [0.5065, 0.35, 0.8939, 0.0]
        upper_right= [1.1, 1.1, 0.0, 0.0]
        lower_left = [1.1, 1.1, 0.8939, -0.8938]
        lower_right= [0.5065, 0.35, 0.0, -0.8939]       
       
        QuadrantTestCase.__init__(self, "Case4", upper_left, upper_right, lower_left, lower_right)

     
class ShearLayer(TestCase):
    
    def __init__(self):
    
        boundary_conditions = [BoundaryPeriodic("left"), BoundaryPeriodic("right"), BoundaryPeriodic("top"), BoundaryPeriodic("bottom")]
    
        TestCase.__init__(self,"ShearLayer", boundary_conditions)
        
    def get_initial_condition(self, grid):
        
        X = grid.X
        Y = grid.Y
        nx = grid.nx
        ny = grid.ny
        ngc = grid.ngc
        Wtemp = np.zeros((nx, ny, 4))

        #width of the stripe
        stripe_w = np.amax(Y) * 0.02
        #center of the stripe
        yc = np.amax(Y) * 0.5


        for i in range(nx):
            for j in range(ny):
                x = X[i,j]
                y = Y[i,j]
                #add a small random component to uy to speedup the breakup
                random = np.random.rand(1)[0] * 0.01
                
                if (abs(y-yc) < stripe_w):
                    
                    Wtemp[i,j, :] = [1.0, 2.5, -0.5, random]  
                else:
                    
                    Wtemp[i,j, :] = [1.0, 2.5, 0.5, random] 

        W = np.zeros((nx + 2*ngc, ny + 2*ngc, 4))
        W[ngc:-ngc, ngc:-ngc, :] = Wtemp
        self.apply_neumann_all(W, grid)
        return W
