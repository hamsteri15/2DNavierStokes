import numpy as np
import pylab as pl
import sys
sys.path.append('../differentiators')
sys.path.append('../time_schemes')
sys.path.append('../common')
sys.path.append('../eq_types/euler')
sys.path.append('../boundary_conditions')
sys.path.append('../filters')
from centered import CDN
from pade import Pade
from weno5 import Weno5
from grid_1d import Grid
from shocktubecalc import sod #analytic solution
from euler import Euler



from bogey_bailly import BogeyBailly
from laplacian import Laplacian
from rk4 import Rk4
from rk3 import Rk3
from rk1 import Rk1


from boundary_neumann import BoundaryNeumann
from boundary_periodic import BoundaryPeriodic

"""

W = [rho, p, u]
U = [rho, rho*E, rho*u]
F(U) = [rho*u, u*(rho*E + p), rho*u*u + p]
"""


def main():

    init = "sod_x" #shocktube
    C = 1.0 #grid streching parameter
    end_time = 0.3
    nx = 10
    ngc = 5
    L = 2.0
    gamma = 1.4
    riemann_solver = "VanLeer"
    
    
   
    dt = 0.0001
    adjust_dt = True
    
    
 
    grid = Grid(L, nx, ngc, C)
    
    boundaries = [BoundaryNeumann("left"), BoundaryNeumann("right")]


    #J = jacobian(C, nx)


    #print np.sum(dG)

    #build spatial operator
    #spatial_operator = Pade(grid)
    #spatial_operator = CDN(grid, 6)
    spatial_operator = Weno5(grid)


    #build the spatial filter
    #sigma = 0.001
    #spatial_filter = Laplacian(grid, sigma, 6)

    sigma = 1
    spatial_filter = BogeyBailly(grid, sigma)

    #build time scheme
    time_scheme = Rk1(dt, adjust_dt)
    #time_scheme = Rk3(dt, adjust_dt)

    #build the equation set
    Equation = Euler(gamma, time_scheme, spatial_operator, riemann_solver)


    W = initial_condition(grid, init)
    
    
    iters = 0
    time = 0.0
    while (time < end_time):

        

        for boundary in boundaries:
            boundary.apply(W, grid)

        W = Equation.take_step(W)
        print W[:,2]
        sys.exit()
        if ( (iters % 20 == 0) and (iters != 0)):
            compare_sod_solutions(W, gamma, time, grid)

        iters +=1
        time += time_scheme.dt
        print "time = {0}   dt = {1} ".format(time, time_scheme.dt)


    #plot_primitive(W, time, nx, ngc, x, L, gamma) 
    compare_sod_solutions(W, gamma, time, grid)
    



def plot_primitive(W, time, nx, ngc, x, L, gamma):

    rho = W[ngc:-ngc, 0]
    p = W[ngc:-ngc,1]
    u = W[ngc:-ngc,2]
    
    pl.plot(x, rho, color="black")
    pl.plot(x, u, color="red")
    pl.plot(x, p, color="blue")
    
    pl.ylim(-0.1,1.1)
    pl.xlim(-1,1)
    pl.show()


def compare_sod_solutions(W, gamma, time, grid):

    L = grid.Lx
    nx = grid.nx
    ngc = grid.ngc

    rho, u, p, E = get_analytic_sod(L, nx, gamma, time)
    
    x_analytic = np.linspace(0, L, nx)
    
    pl.plot(x_analytic, rho, color="black", ls="--")
    pl.plot(x_analytic, u, color="red", ls="--")
    pl.plot(x_analytic, p, color="blue", ls="--")
    
    
    rho = W[ngc:-ngc, 0]
    p = W[ngc:-ngc,1]
    u = W[ngc:-ngc,2]
    
    pl.scatter(grid.X, rho, color="black", label="rho")
    pl.scatter(grid.X, u, color="red", label ="u")
    pl.scatter(grid.X, p, color="blue", label="p")
    
    pl.ylim(-0.1,1.1)
    pl.xlim(0,L)
    pl.legend(loc="best")
    pl.show()

    

    


def Neumann_bc(array, ngc):
	
	#left
    for i in range(0, ngc):
        array[i] = array[ngc]
	
	#right
    for i in range(0,ngc):
        array[-i-1] = array[-ngc-1]	
        


def apply_periodic_boundary(f, ngc):

    f_new = f
    f_new[0:ngc] = f[-2*ngc:-ngc]
    f_new[-ngc:] = f[ngc:2*ngc]
    
    return f_new


def initial_condition(grid, init):
    
    X = grid.X
    nx = grid.nx
    ngc = grid.ngc

    Wtemp = np.zeros((nx, 3))
    if (init == "sod_x"):
        left_state=[1.0, 1.0, 0.0]  #rho p u v
        right_state=[0.125, 0.1, 0.0]
        Wtemp[0:int(nx/2), :] = left_state    
        Wtemp[int(nx/2):, :] = right_state
        
    


    W = np.zeros((nx + 2*ngc, 3))
    W[ngc:-ngc, :] = Wtemp


    return W
    







def get_analytic_sod(L, nx, gamma, time=0.):
	
	#Initializes primitive variables for Sod's shock tube problem
	
	
	positions, regions, values = sod.solve(left_state=(1, 1, 0), right_state=(0.1, 0.125, 0.), geometry=(-0.5*L, 0.5*L, 0.), t=time, gamma=gamma, npts=nx)
	
	
	rho0 = values["rho"]	
	u0 = values["u"]
	p0 = values["p"]
	
	rho0[-1] = rho0[-2]
	u0[-1] = u0[-2]
	p0[-1] = p0[-2]
	
	E0 = p0/((gamma-1)*rho0) + 0.5*u0**2 #total energy

	return rho0,u0,p0,E0
































main()
