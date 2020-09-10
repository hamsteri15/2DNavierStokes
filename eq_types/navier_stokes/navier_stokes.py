import numpy as np
import pylab as pl
import sys

from euler import Euler

#Extends the Euler class with diffusion scheme
class NavierStokes(Euler):

    def __init__(self, grid_parameters, time_parameters, solver_parameters, fluid_properties):
        
        Euler.__init__(self, grid_parameters, time_parameters, solver_parameters, fluid_properties)
        
        self.spatial_operator_diffusion = solver_parameters.get_diffusion_scheme(self.grid)
        self.viscosity = solver_parameters.get_viscosity_method(fluid_properties)
        
    def dU(self, U):

        return self.dU_convection_2d(U) + self.dU_diffusion_2d(U)



    def tauxx(self, mu, du_dx, dv_dy, lamb):

        ret = 2.0 * mu * du_dx + lamb * (du_dx + dv_dy)
        return ret

    def tauyy(self, mu, du_dx, dv_dy, lamb):

        ret = 2.0 * mu * dv_dy + lamb * (du_dx + dv_dy)
        return ret

    def tauxy(self, mu, du_dy, dv_dx):
        
        ret = mu * (du_dy + dv_dx)
        return ret

    def dtauxx_dx(self, mu, dmu_dx, du_dx, dv_dy, du_dxdx, dv_dydx, lamb, dlamb_dx):
        
        ret = 2.0 * dmu_dx * du_dx + 2.0 * mu * du_dxdx + dlamb_dx*(du_dx + dv_dy)
        ret += lamb * (du_dxdx + dv_dydx) 
        return ret

    def dtauyy_dy(self, mu, dmu_dy, dv_dy, dv_dydy, du_dx, du_dxdy, lamb, dlamb_dy):

        ret = 2.0 * dmu_dy * dv_dy + 2.0 * mu * dv_dydy + dlamb_dy * (du_dx + dv_dy)
        ret += lamb * (du_dxdy + dv_dydy)      
        return ret

    def dtauxy_dx(self, mu, dmu_dx, du_dy, dv_dx, du_dydx, dv_dxdx):

        ret = dmu_dx * (du_dy + dv_dx) + mu * (du_dydx + dv_dxdx)
        return ret

    def dtauyx_dy(self, mu, dmu_dy, du_dy, dv_dx, du_dydy, dv_dxdy):

        ret = dmu_dy * (du_dy + dv_dx) + mu * (du_dydy + dv_dxdy)
        return ret


    def dU_diffusion_2d(self, U):


        d_dx = self.spatial_operator_diffusion.d_dx
        d_dy = self.spatial_operator_diffusion.d_dx
        d_dxdx = self.spatial_operator_diffusion.d_dxdx
        d_dydy = self.spatial_operator_diffusion.d_dydy

        #this step is not necessary since primitives are already available in take_step()
        W = self.conservative_to_primitive(U)

        #viscous residual
        residual = np.zeros(U.shape)

        rho = W[:,:,0]
        p = W[:,:,1]
        u = W[:,:,2]
        v = W[:,:,3]
        T = self.viscosity.temperature(p, rho)

        mu = self.viscosity.mu(T)

        #print np.mean(mu)
        #print np.mean(T)
        

        dmu_dx = self.viscosity.dmu_dx(self.spatial_operator_diffusion, mu)
        dmu_dy = self.viscosity.dmu_dy(self.spatial_operator_diffusion, mu)
        lamb = -(2./3.) * mu
        dlamb_dx = -(2./3.) * dmu_dx
        dlamb_dy = -(2./3.) * dmu_dy
        
        """
        X = self.spatial_operator_diffusion.grid.X
        Y = self.spatial_operator_diffusion.grid.Y
        ngc = self.spatial_operator_diffusion.grid.ngc

        pl.pcolor(X,Y,dmu_dx[ngc:-ngc, ngc:-ngc])
        pl.colorbar()
        pl.show()
        """
        
        du_dx = d_dx(u)
        du_dy = d_dy(u)
        du_dxdx = d_dxdx(u)
        du_dydy = d_dydy(u)
        du_dxdy = d_dy(d_dx(u))
        du_dydx = du_dxdy #analytically ok, discretely not sure...


        dv_dy = d_dy(v)
        dv_dx = d_dx(v)
        dv_dydy = d_dydy(v)
        dv_dxdx = d_dxdx(v)
        dv_dydx = d_dx(d_dy(v))
        dv_dxdy = dv_dydx #analytically ok, discretely not sure...

        
        drho_dx = d_dx(rho)
        dp_dx = d_dx(p)
        drho_dxdx = d_dxdx(rho)
        dp_dxdx = d_dxdx(p)
        


        drho_dy = d_dy(rho)
        dp_dy = d_dy(p)
        drho_dydy = d_dydy(rho)
        dp_dydy = d_dydy(p)




        Dtauxx_dx = self.dtauxx_dx(mu, dmu_dx, du_dx, dv_dy, du_dxdx, dv_dydx, lamb, dlamb_dx)
        Dtauyy_dy = self.dtauyy_dy(mu, dmu_dy, dv_dy, dv_dydy, du_dx, du_dxdy, lamb, dlamb_dy)

        Dtauxy_dx = self.dtauxy_dx(mu, dmu_dx, du_dy, dv_dx, du_dydx, dv_dxdx)
        Dtauyx_dy = self.dtauyx_dy(mu, dmu_dy, du_dy, dv_dx, du_dydy, dv_dxdy)


        tauxx = self.tauxx(mu, du_dx, dv_dy, lamb)
        tauyy = self.tauxx(mu, du_dx, dv_dy, lamb)
        tauxy = self.tauxy(mu, du_dy, dv_dx)
        tauyx = tauxy #this is correct for sure



        Dtauvx_dx = u * Dtauxx_dx + v * Dtauxy_dx + tauxx * du_dx + tauxy * dv_dy
        Dtauvy_dy = u * Dtauyx_dy + v * Dtauyy_dy + tauyy * dv_dy + tauyx * du_dy

        one_rho = 1./rho
        one_rho2 = 1./(rho*rho)
        one_rho3 = 1./(rho*rho*rho)


        #derivative of heat flux x-dir
        temp1 = dmu_dx * (one_rho * dp_dx - p * one_rho2 * drho_dx) 
        temp2 =  mu * (2 * one_rho2 * drho_dx * dp_dx + one_rho * dp_dxdx + one_rho3 * p * drho_dx * drho_dx - one_rho2 * p * drho_dxdx)
        temp3 = self.gamma / (self.Pr * (self.gamma - 1.0))

        Dqx_dx = -temp3 * (temp1 + temp2)

        #derivative of heat flux y-dir
        temp1 = dmu_dy * (one_rho * dp_dy - p * one_rho2 * drho_dy)
        temp2 = mu * (2 * one_rho2 * drho_dy * dp_dy + one_rho * dp_dydy + one_rho3 * p * drho_dy * drho_dy - one_rho2 * p * drho_dydy)
        Dqy_dy = -temp3 * (temp1 + temp2)

        x_mom = -Dtauxx_dx -Dtauyx_dy
        y_mom = -Dtauyy_dy -Dtauxy_dx
        ene =   -Dtauvx_dx - Dtauvy_dy + Dqx_dx + Dqy_dy  #+heat fluxes

        #convection residual is positive!
        residual[:,:,1] = ene
        residual[:,:,2] = x_mom
        residual[:,:,3] = y_mom
        

        """
        #x_mom = (4./3.) * d_dxdx(u) + d_dydy(u) + (1./3.) * d_dx(d_dy(v))
        #y_mom = (4./3.) * d_dydy(v) + d_dxdx(v) + (1./3.) * d_dx(d_dy(u))
        x_mom = (4./3.) * du_dxdx + du_dydy + (1./3.) * dv_dxdy
        y_mom = (4./3.) * dv_dydy + dv_dxdx + (1./3.) * du_dydx

        residual[:,:,2] = self.viscosity.mu() * x_mom
        residual[:,:,3] = self.viscosity.mu() * y_mom
        """

        return residual








        
       
    
