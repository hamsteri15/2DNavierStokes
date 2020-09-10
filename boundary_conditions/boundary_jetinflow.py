import numpy as np
import sys



class BoundaryJetinflow:

    def __init__(self, location, diameter, ux):
        
        #either "left, right, top or bottom"
        self.location = location
        self.D = diameter
        self.ux = ux


        #only to left boundary
        if (location == "left"):
            self.apply = self.apply_left
          

        else:
            print "invalid boundary location: " + location
            sys.exit()

    def apply_left(self, W, grid):
        
        ngc = grid.ngc
        if (grid.is_1d()):
            print "no jet inflow boundary for 1d-meshes!"
            sys.exit()

        else:

            ngc = grid.ngc
            ny = grid.ny
            yc = 0.5*grid.Ly
            #left (this is the inlet)
            for j in range(ny):
                for i in range(0, ngc):
                    jj = j + ngc
                    y = grid.Y[0,j]

                    if (abs(y - yc) < 0.5*self.D):
                        W[i, jj, :] = [1.0, 1.0, self.ux, 0.0]
                    else:
                        W[i, jj, :] = [1.0, 1.0, 0.0, 0.0]
                  

