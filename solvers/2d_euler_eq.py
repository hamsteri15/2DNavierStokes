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
    
    
    grid_p = GridParameters(nx=100, ny=100, ngc=3, Cx=0.0, Cy=0.0, Lx=1.0, Ly=1.0)
    time_p = TimeParameters(dt = 0.0001, end_t = 0.3, time_scheme_name="RK-3", adjust_dt = True)
    
    euler_p = EulerSolverParameters(convection_scheme_name="Weno", \
                                    convection_flux_type="LaxFriedrichs", \
                                    filter_name="None",
                                    filter_sigma = 0.0)
                                    
                                    
                                    
                                    
    fluid_p = FluidPropertiesAir()
    
    case = Case3()
    
    Equation = Euler(grid_p, time_p, euler_p, fluid_p)
    W = case.get_initial_condition(Equation.grid)
    

    iters = 0
    time = 0.0
    while (time < time_p.end_t):

        case.apply_boundary_conditions(W, Equation.grid)
 
        W = Equation.take_step(W)
        
        
        if ( (iters % 10 == 0) and (iters != 0)):
            #pass
            
            plot_density(Equation.grid, W)
            #plot_pressure(grid, W)
            
        time += Equation.dt
        iters+=1
        print "time = {0}   dt = {1} ".format(time, Equation.dt)
    
    
    
    
    print "done"





    
def plot_pressure(grid, W):
    
    ngc = grid.ngc
    X = grid.X
    Y = grid.Y

    pl.pcolor(X,Y,W[ngc:-ngc, ngc:-ngc, 1])
    pl.colorbar()
    pl.xlim(np.amin(X), np.amax(X))
    pl.ylim(np.amin(Y), np.amax(Y))
    pl.show()



def plot_density(grid, W):
    
    ngc = grid.ngc
    X = grid.X
    Y = grid.Y

    pl.pcolor(X,Y,W[ngc:-ngc, ngc:-ngc, 0], edgecolors="black" )
    #pl.pcolor(X,Y,W[ngc:-ngc, ngc:-ngc, 0] )
    pl.colorbar()
    pl.xlim(np.amin(X), np.amax(X))
    pl.ylim(np.amin(Y), np.amax(Y))
    pl.show()

def compare_sod_solutions(W, time, nx, ngc, x, L, gamma):

    rho, u, p, E = get_analytic_sod(L, nx, gamma, time)
    
    x_analytic = np.linspace(-0.5*L, 0.5*L, nx)
    
    pl.plot(x_analytic, rho, color="black", ls="--")
    pl.plot(x_analytic, u, color="red", ls="--")
    pl.plot(x_analytic, p, color="blue", ls="--")
    
    
    rho = W[ngc:-ngc, 0]
    p = W[ngc:-ngc,1]
    u = W[ngc:-ngc,2]
    
    pl.plot(x, rho, color="black")
    pl.plot(x, u, color="red")
    pl.plot(x, p, color="blue")
    
    pl.ylim(0,1)
    pl.xlim(-1,1)
    pl.show()




def Neumann_bc(array, ngc):
	
	#left
    for i in range(0, ngc):
        array[i, :, :] = array[ngc, :, :]
	
	#right
    for i in range(0, ngc):
        array[-i-1, :, :] = array[-ngc-1, :, :]	
        
    
    #top
    for j in range(0, ngc):
        array[:, j, :] = array[:, ngc, :]
	
	#bottom
    for j in range(0, ngc):
        array[:, -j-1, :] = array[:, -ngc-1, :]	


def periodic_bc(array, ngc, nx, ny):

    #left
    for i in range(0, ngc):
        for j in range(ngc, ny + ngc):
            for q in range(4):
                array[i, j, q] =  array[nx + i, j, q]
	
    
	#right
    for i in range(0, ngc):
        for j in range(ngc, ny + ngc):
            for q in range(4):
                array[nx + ngc + i, j, q] = array[i + ngc, j, q]	
        
    
    #top
    for i in range(ngc, nx + ngc):
        for j in range(0, ngc):
            for q in range(4):
                array[i, j, q] = array[i, ny + j, q]
	
	#bottom
    for i in range(ngc, nx + ngc):
        for j in range(0, ngc):
            for q in range(4):
                array[i, ny + ngc + j, q] = array[i, j + ngc, q]
        	
    #print np.round(array[:, :, 0],2)
    #print "---"
    #print array[ngc + 1, :, 0]
    #sys.exit()

def initial_condition(grid, init):
    
    X = grid.X
    Y = grid.Y
    nx = grid.nx
    ny = grid.ny
    ngc = grid.ngc

    Wtemp = np.zeros((nx, ny, 4))
    if (init == "sod_x"):
        left_state=[1.0, 1.0, 0.0, 0.0]  #rho p u v
        right_state=[0.125, 0.1, 0.0, 0.0]
        Wtemp[0:int(nx/2), :, :] = left_state    
        Wtemp[int(nx/2):, :, :] = right_state
        W = np.zeros((nx + 2*ngc, ny + 2*ngc, 4))
        W[ngc:-ngc, ngc:-ngc, :] = Wtemp  
        return W
    if (init == "sod_y"):
        left_state=[1.0, 1.0, 0.0, 0.0]  #rho p u v
        right_state=[0.125, 0.1, 0.0, 0.0]
        Wtemp[:, 0:int(ny/2), :] = left_state    
        Wtemp[:, int(ny/2):, :] =  right_state
        W = np.zeros((nx + 2*ngc, ny + 2*ngc, 4))
        W[ngc:-ngc, ngc:-ngc, :] = Wtemp  
        return W
    
    if (init == "case3"):
        upper_left=  [0.5323, 0.3, 1.206, 0.0]  #rho p u v
        upper_right= [1.5, 1.5, 0.0, 0.0]
        lower_left = [0.138, 0.029, 1.206, -1.206]
        lower_right= [0.5323, 0.3, 0.0, -1.206]
        
    if (init == "case4"):
        upper_left=  [0.5065, 0.35, 0.8939, 0.0]  #rho p u v
        upper_right= [1.1, 1.1, 0.0, 0.0]
        lower_left = [1.1, 1.1, 0.8939, -0.8938]
        lower_right= [0.5065, 0.35, 0.0, -0.8939]


    if (init == "shear_layer"):

        return set_shear_layer(X, Y, nx, ny, ngc)       


        #std::vector<double> NW = {0.5323, 0.3, 1.206, 0.0, 0.0};
	    #std::vector<double> NE = {1.5, 1.5, 0.0, 0.0, 0.0};
	    #std::vector<double> SW = {0.138, 0.029, 1.206, -1.206, 0.0};
	    #std::vector<double> SE = {0.5323, 0.3, 0.0, -1.206, 0.0};

        
    Wtemp[0:int(nx/2), 0:int(ny/2), :] = upper_left
    Wtemp[int(nx/2):, 0:int(ny/2), :] = upper_right
    Wtemp[0:int(nx/2):, int(ny/2):, :] = lower_left    
    Wtemp[int(nx/2):, int(ny/2):, :] =  lower_right
        


    W = np.zeros((nx + 2*ngc, ny + 2*ngc, 4))
    W[ngc:-ngc, ngc:-ngc, :] = Wtemp


    return W
    


def set_shear_layer(X, Y, nx, ny, ngc):

    Wtemp = np.zeros((nx, ny, 4))
    
    stripe_w = 0.25

    for i in range(nx):
        for j in range(ny):
            x = X[i,j]
            y = Y[i,j]

            #no fucking idea why the velocity components have to be in wrong order
            if (abs(y) < stripe_w):
                #Wtemp[i,j] = [2.0, 2.5, -0.5, 0.0*np.sin(np.pi*x)] #0.01
                Wtemp[i,j, :] = [2.0, 2.5, 0.01*np.sin(np.pi*x), -0.5] #0.01
            else:
                #Wtemp[i,j] = [1.0, 2.5, 0.5, 0.0*np.sin(np.pi*x)]
                Wtemp[i,j, :] = [1.0, 2.5, 0.01*np.sin(np.pi*x), 0.5] #0.01

    W = np.zeros((nx + 2*ngc, ny + 2*ngc, 4))
    W[ngc:-ngc, ngc:-ngc, :] = Wtemp
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
