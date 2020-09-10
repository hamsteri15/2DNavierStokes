import numpy as np
import sys
import copy
sys.path.append('../differentiators')
from centered import CDN


class Laplacian:

    def __init__(self, grid, sigma, order=6):
        
        self.filter_width = sigma

       
        self.my_grid = grid
        self.spatial_operator = CDN(self.my_grid, order)
        

    
        



    def filter(self, u):
        
    

        if (self.my_grid.is_1d()):
            laplacian = self.spatial_operator.d_dxdx(u)
            multiplier = self.my_grid.dx * self.my_grid.dx / (np.pi*np.pi)
            u_new = u + self.filter_width * multiplier * laplacian

        if(self.my_grid.is_2d()):
            
            multiplier_x = self.filter_width * self.my_grid.dx * self.my_grid.dx / (np.pi*np.pi)
            multiplier_y = self.filter_width * self.my_grid.dy * self.my_grid.dy / (np.pi*np.pi)

            u_filtered = multiplier_x * self.spatial_operator.d_dxdx(u) + multiplier_y * self.spatial_operator.d_dydy(u)
            u_new = u + u_filtered

        
        
        

        return u_new
        

