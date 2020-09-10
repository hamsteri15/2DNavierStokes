import pylab as pl



def plot_v(grid, W):
    
    ngc = grid.ngc
    X = grid.X
    Y = grid.Y

    pl.pcolor(X,Y,W[ngc:-ngc, ngc:-ngc, 3])
    pl.colorbar()
    pl.xlim(np.amin(X), np.amax(X))
    pl.ylim(np.amin(Y), np.amax(Y))
    pl.show()

def plot_u(grid, W):
    
    ngc = grid.ngc
    X = grid.X
    Y = grid.Y

    pl.pcolor(X,Y,W[ngc:-ngc, ngc:-ngc, 2])
    pl.colorbar()
    pl.xlim(np.amin(X), np.amax(X))
    pl.ylim(np.amin(Y), np.amax(Y))
    pl.axes().set_aspect('equal', 'datalim')
    #pl.show()


def plot_pressure(grid, W):
    
    ngc = grid.ngc
    X = grid.X
    Y = grid.Y

    pl.pcolor(X,Y,W[ngc:-ngc, ngc:-ngc, 1])
    pl.colorbar()
    pl.xlim(np.amin(X), np.amax(X))
    pl.ylim(np.amin(Y), np.amax(Y))
    #pl.show()



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


def plot_density_grad(grid, W):
    
    ngc = grid.ngc
    X = grid.X
    Y = grid.Y


    diff = CDN(grid, 6)

    rho = W[:,:,0]
    v = W[:,:,3]
    drho_dx  = diff.d_dx(rho)
    drho_dy  = diff.d_dy(rho)
    grad = np.sqrt(drho_dx * drho_dx + drho_dy * drho_dy)


    pl.pcolor(X,Y,grad[ngc:-ngc, ngc:-ngc], edgecolors="black")
    pl.colorbar()
    pl.xlim(np.amin(X), np.amax(X))
    pl.ylim(np.amin(Y), np.amax(Y))
    pl.axes().set_aspect('equal', 'datalim')
    #pl.show()

def plot_vorticity(grid, W, showGrid=False, showImage=False):
    
    import sys
    import numpy as np
    sys.path.append('../differentiators')
    from centered import CDN
    
    ngc = grid.ngc
    X = grid.X
    Y = grid.Y


    diff = CDN(grid, 2)

    u = W[:,:,2]
    v = W[:,:,3]
    du_dy  = diff.d_dy(u)
    dv_dx  = diff.d_dx(v)
    omega = du_dy - dv_dx

    if (showGrid):
        pl.pcolor(X,Y,omega[ngc:-ngc, ngc:-ngc], edgecolors="black")
    else:
        pl.pcolor(X,Y,omega[ngc:-ngc, ngc:-ngc])
    pl.colorbar()
    pl.xlim(np.amin(X), np.amax(X))
    pl.ylim(np.amin(Y), np.amax(Y))
    if (showImage):
        pl.show()
