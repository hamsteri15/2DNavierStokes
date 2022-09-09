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
 
    W = Equation.take_step(W)


    
    
    
    print ("done")




























main()
