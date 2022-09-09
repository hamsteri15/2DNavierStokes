import sys

sys.path.append('../differentiators')
sys.path.append('../time_schemes')
sys.path.append('../eq_types/navier_stokes')
sys.path.append('../eq_types/euler')
sys.path.append('../filters')

from centered import CDN
from pade import Pade
from weno5 import Weno5
from sutherland_viscosity import SutherlandViscosity
from constant_viscosity import ConstantViscosity
from bogey_bailly import BogeyBailly
from laplacian import Laplacian
from filter_none import FilterNone
from rk4 import Rk4
from rk3 import Rk3
from rk1 import Rk1
from laxfriedrichs_flux import LaxFriedrichs
from vanleer_flux import VanLeer
from physical_flux import PhysicalFlux


class GridParameters:
    def __init__(self, nx, ny, ngc, Cx, Cy, Lx, Ly):
        self.nx = nx
        self.ny = ny
        self.ngc = ngc
        self.Cx = Cx
        self.Cy = Cy
        self.Lx = Lx
        self.Ly = Ly


class TimeParameters:

    def __init__(self, dt, end_t, time_scheme_name, adjust_dt = True):
    
        
        self.dt = dt
        self.end_t = end_t
        self.time_scheme_name = time_scheme_name
        self.adjust_dt = adjust_dt
        
        
    def get_time_scheme(self):
    
        if (self.time_scheme_name == "RK-1"):
    
            return Rk1()
            
        elif (self.time_scheme_name == "RK-3"):
    
            return Rk3()
        
        elif (self.time_scheme_name == "RK-4"):
    
            return Rk4()
            
        else:
            print ("Invalid time scheme")
            sys.exit()
        
class EulerSolverParameters:
    
    def __init__(self, convection_scheme_name, convection_flux_type, filter_name, filter_sigma):
        
        
        self.convection_scheme_name = convection_scheme_name
        self.convection_flux_type = convection_flux_type
        self.filter_name = filter_name
        self.filter_sigma = filter_sigma
        
    def get_convection_scheme(self, grid):
        if (self.convection_scheme_name == "Weno"):
            return Weno5(grid)
            
        elif (self.convection_scheme_name == "CD-2"):
            
            return CDN(grid, 2)
            
        elif (self.convection_scheme_name == "CD-4"):
            
            return CDN(grid, 4)
            
        elif (self.convection_scheme_name == "CD-6"):
            
            return CDN(grid, 6)
            
        elif (self.convection_scheme_name == "Pade"):
            
            return Pade(grid)
            
        else:
            print ("Invalid convection scheme")
            sys.exit()
            
    def get_convection_flux(self, gamma):
    
        if (self.convection_flux_type == "Physical"):
            return PhysicalFlux(gamma)
    
        elif (self.convection_flux_type == "VanLeer"):
            return VanLeer(gamma)
    
        elif (self.convection_flux_type == "LaxFriedrichs"):
            return LaxFriedrichs(gamma)
            
        else:
            print ("Invalid flux type")
            sys.exit()
           
            
    def get_filter(self, grid):
        
        if (self.filter_name == "Bogey-Bailly"):
            return BogeyBailly(grid, self.filter_sigma)
            
        elif (self.filter_name == "Laplacian"):
            return Laplacian(grid, self.filter_sigma)
            
        else:
            return FilterNone()
            

class NavierStokesSolverParameters(EulerSolverParameters):
    
    
    def __init__(self, convection_scheme_name, convection_flux_type, filter_name, filter_sigma, diffusion_scheme_name, viscosity_type):
        
        EulerSolverParameters.__init__(self,convection_scheme_name, convection_flux_type, filter_name, filter_sigma)
        
        self.diffusion_scheme_name = diffusion_scheme_name
        self.viscosity_type = viscosity_type
        
    def get_diffusion_scheme(self, grid):
        
        
        if (self.diffusion_scheme_name == "Pade"):
                return Pade(grid)
                
        elif (self.diffusion_scheme_name == "CD-2"):
            
            return CDN(grid, 2)
            
        elif (self.diffusion_scheme_name == "CD-4"):
            
            return CDN(grid, 4)
            
        elif (self.diffusion_scheme_name == "CD-6"):
            
            return CDN(grid, 6)
            
        else:
            print ("Invalid convection scheme")
            sys.exit()


    def get_viscosity_method(self, fluid):
        if (self.viscosity_type == "Constant"):
            return ConstantViscosity(fluid.mu0, fluid.R) 
        
        elif (self.viscosity_type == "Sutherland"):
            return SutherlandViscosity(fluid.mu0, fluid.T0, fluid.C, fluid.R)
        
        else:
            print ("Invalid viscosity type")
            sys.exit()
            
            
            
            
            
            
            
            
            
            
