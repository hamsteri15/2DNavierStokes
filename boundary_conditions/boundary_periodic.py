import numpy as np
import sys



class BoundaryPeriodic:

    def __init__(self, location):
        
        #either "left, right, top or bottom"
        self.location = location

        if (location == "left"):
            self.apply = self.apply_left
        elif (location == "right"):
            self.apply = self.apply_right
        elif (location == "top"):
            self.apply = self.apply_top
        elif (location == "bottom"):
            self.apply = self.apply_bottom    

        else:
            print "invalid boundary location: " + location
            sys.exit()

    def apply_left(self, W, grid):
        
        ngc = grid.ngc
        nx = grid.nx
        ny = grid.ny
        if (grid.is_1d()):
            for i in range(0, ngc):
                W[i, :] =  W[nx + i, :] 

        else:
            
            for i in range(0, ngc):
                for j in range(ngc, ny + ngc):
                        W[i, j, :] =  W[nx + i, j, :]        



    def apply_right(self, W, grid):
        
        ngc = grid.ngc
        nx = grid.nx
        ny = grid.ny
        if (grid.is_1d()):
            for i in range(0, ngc):
                W[nx + ngc + i, :] = W[i + ngc, :] 

        else:
            
            for i in range(0, ngc):
                for j in range(ngc, ny + ngc):
                        W[nx + ngc + i, j, :] = W[i + ngc, j, :]      
        

    def apply_top(self, W, grid):
        
        ngc = grid.ngc
        nx = grid.nx
        ny = grid.ny
        if (grid.is_1d()):
            print "applying top boundary for 1d grid, DONT!"
            sys.exit()	 

        else:
            for i in range(ngc, nx + ngc):
                for j in range(0, ngc):
                        W[i, j, :] = W[i, ny + j, :]	
        

    def apply_bottom(self, W, grid):
        
        ngc = grid.ngc
        nx = grid.nx
        ny = grid.ny
        if (grid.is_1d()):
            print "applying bottom boundary for 1d grid, DONT!"
            sys.exit()	
            
        else:
            #bottom
            for i in range(ngc, nx + ngc):
                for j in range(0, ngc):
                        W[i, ny + ngc + j, :] = W[i, j + ngc, :]



    