import numpy as np
import sys



class BoundaryNeumann:

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
            print ("invalid boundary location: " + location)
            sys.exit()

    def apply_left(self, W, grid):
        
        ngc = grid.ngc
        if (grid.is_1d()):
            for i in range(0, ngc):
                W[i, :] = W[ngc, :]

        else:
            #left
            for i in range(0, ngc):
                W[i, :, :] = W[ngc, :, :]        



    def apply_right(self, W, grid):
        
        ngc = grid.ngc

        if (grid.is_1d()):
            for i in range(0, ngc):
                W[-i-1, :] = W[-ngc-1, :]	
            
        else:
            #right
            for i in range(0, ngc):
                W[-i-1, :, :] = W[-ngc-1, :, :]	     
        

    def apply_top(self, W, grid):
        
        ngc = grid.ngc

        if (grid.is_1d()):
            print ("applying top boundary for 1d grid, DONT!")
            sys.exit()	
            
        else:
            #top
            for j in range(0, ngc):
                W[:, j, :] = W[:, ngc, :]	
        

    def apply_bottom(self, W, grid):
        
        ngc = grid.ngc

        if (grid.is_1d()):
            print ("applying bottom boundary for 1d grid, DONT!")
            sys.exit()	
            
        else:
            #bottom
            for j in range(0, ngc):
                W[:, -j-1, :] = W[:, -ngc-1, :]	



    
