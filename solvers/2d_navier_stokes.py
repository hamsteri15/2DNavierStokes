import numpy as np
import pylab as pl
import sys

sys.path.append('../common')
sys.path.append('../eq_types/navier_stokes')
sys.path.append('../eq_types/euler')
sys.path.append('../boundary_conditions')
sys.path.append('../filters')
sys.path.append('../test_cases')

from shocktubecalc import sod #analytic solution
from navier_stokes import NavierStokes
from solver_parameters import *
from fluid_properties import *
from test_cases import *
from plot import *

def main():

    #Cs = np.linspace(0.01, 3., 4)
    
    
    grid_p = GridParameters(nx=250, ny=250, ngc=3, Cx=1.0, Cy=3.0, Lx=1.0, Ly=1.0)
    time_p = TimeParameters(dt = 0.0001, end_t = 2, time_scheme_name="RK-3", adjust_dt = True)
    
    
    ns_p = NavierStokesSolverParameters(convection_scheme_name="CD-6", \
                                        convection_flux_type="Physical", \
                                        filter_name="Laplacian",
                                        filter_sigma = 0.05, \
                                        diffusion_scheme_name="CD-6", \
                                        viscosity_type = "Constant")
                                
    
    fluid_p = FluidPropertiesAir()
    
    case = ShearLayer()
    
    Equation = NavierStokes(grid_p, time_p, ns_p, fluid_p)
    W = case.get_initial_condition(Equation.grid)
    

    iters = 0
    time = 0.0
    while (time < time_p.end_t):

        case.apply_boundary_conditions(W, Equation.grid)
 
        W = Equation.take_step(W)
        
        
        if ( (iters % 100 == 0) and (iters != 0)):
            #pass
            
            plot_vorticity(Equation.grid, W)
            pl.savefig("vorticity_{0}.png".format(iters))
            #pl.show()
            pl.clf()
            
        time += Equation.dt
        iters+=1
        print "time = {0}   dt = {1} ".format(time, Equation.dt)
    
 
    
    
    print "done"



































main()
