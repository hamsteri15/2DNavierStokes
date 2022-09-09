import numpy as np
import pylab as pl
import sys

sys.path.append('../common')
sys.path.append('../eq_types/euler')
sys.path.append('../boundary_conditions')
sys.path.append('../filters')
sys.path.append('../test_cases')

from shocktubecalc import sod #analytic solution
from euler import Euler
from solver_parameters import *
from fluid_properties import *
from test_cases import *
from grid_2d import Grid

def main():

    #Cs = np.linspace(0.01, 3., 4)
    
    
    grid_p = GridParameters(nx=5, ny=1, ngc=3, Cx=0.0, Cy=0.0, Lx=1.0, Ly=1.0)
    time_p = TimeParameters(dt = 0.0001, end_t = 0.3, time_scheme_name="RK-1", adjust_dt = False)
    
    euler_p = EulerSolverParameters(convection_scheme_name="Weno", \
                                    convection_flux_type="LaxFriedrichs", \
                                    filter_name="None",
                                    filter_sigma = 0.0)
                                    
                                    
             
    fluid_p = FluidPropertiesAir()
    
    case = Sod_x("sod_x")
    
    Equation = Euler(grid_p, time_p, euler_p, fluid_p)
    W = case.get_initial_condition(Equation.grid)
    
    case.apply_boundary_conditions(W, Equation.grid)
 
    eig_x, eig_y = Equation.get_eigenvalues(W)

    Equation.eig_x = eig_x
    Equation.eig_y = eig_y

    U = Equation.primitive_to_conservative(W)

    """
    flux_method = PhysicalFlux(Equation.gamma)

    Fx = flux_method.flux(U, "x")
    print (Fx)

    """

    """
    fl, fr = Equation.flux_method.flux(U, Equation.eig_x, direction="x")

    print (fl)
    #print (fr[:,:,2])
   

    """

    U = Equation.primitive_to_conservative(W)
    dU = Equation.dU_convection_2d(U)
    
    rho = dU[:, 3, 0]
    rhoE = dU[:, 3, 1]
    rhou = dU[:, 3, 2]
    rhov = dU[:, 3, 3]
    
    print ("rho: ")
    print (rho)

    """
    print ("rhoE: ")
    print (rhoE)

    print ("rhou: ")
    print (rhou)

    print ("rhov: ")
    print (rhov)
    """
    


    
    
    
    print ("done")




























main()
