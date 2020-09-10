import numpy as np
from numpy.linalg import inv
from scipy.sparse import csr_matrix
import sys #exit


class Pade:
    """
    Generic Pade differentiation class. The nx (and ny) parameters passed to the constructor _should_ include ghost cells,
    otherwise the boundaries will be treated badly.
    """
    

    def __init__(self, grid):
        
        
        
        
        self.grid = grid

        #Pade coeffs for first derivative
        self.alpha_d1 = 1./4.; self.a_d1 =  (2./3.)*(self.alpha_d1 + 2.0);
        #Pade coeffs for the second derivative
        self.alpha_d2 = 1./10.; self.a_d2 = (4./3.)*(1.0 - self.alpha_d2); 



        self.order = 4
        self.name = "Pade-{0}".format(self.order)

        
        

        #1D-Pade operators
        nx = grid.nx + grid.ngc*2 
        self.dx_matrix = self.create_matrix_d(nx)
        self.ddx_matrix = self.create_matrix_dd(nx)


        


        if (self.grid.is_2d()):
            ny = grid.ny + grid.ngc*2
            self.dy_matrix = self.create_matrix_d(ny)
            self.ddy_matrix = self.create_matrix_dd(ny)

        



    def create_matrix_d(self, N):
        """
        matrix for first derivative
        """
        A = np.eye(N,N,k=-1)*self.alpha_d1 + np.eye(N,N)*1 + np.eye(N,N,k=1)*self.alpha_d1        
        B = np.eye(N,N,k=1)*1. + np.eye(N,N,k=-1)*-1. #f_i+1 - f_i-1
        

        

        #A, B, G = self.set_periodic_boundaries(N,A,B,G)
        #A, B, G = self.set_neumann_bc(N,A,B,G)
        #print B
        #sys.exit()


        B = B*(self.a_d1/2.0)
        
        temp = inv(A) 
        rhs = B

        operator = np.dot(temp,rhs)
        return operator


    def create_matrix_dd(self, N):
        """
        matrix for second derivative
        """
        A = np.eye(N,N,k=-1)*self.alpha_d2 + np.eye(N,N)*1 + np.eye(N,N,k=1)*self.alpha_d2        
        B = np.eye(N,N,k=1)*1. + np.eye(N,N)*(-2.0) + np.eye(N,N,k=-1)*1.
        
        
        

        A, B = self.set_periodic_boundaries(N,A,B)
        #A, B, G = self.set_neumann_bc(N,A,B,G)
        #print B
        #sys.exit()


        B = B*self.a_d2
        
        temp = inv(A) 
        rhs = B

        operator = np.dot(temp,rhs)
        return operator



    def set_neumann_bc(self, N, A, B, G):
        """
        
        """
        
        #B[0,0] = 0.5 * (3.0 + self.alfa)#this should be (3 + alpha)/2 according to lele
        #B[N-1,0] = 1.
        B[0,0] = 1.0
        B[0,1] = 0.0
        
        #print B
        #sys.exit()


        Anew = A; Bnew = B; Gnew = G
        
    
    
        return Anew, Bnew, Gnew





    def set_periodic_boundaries(self, N, A, B):
    
        
        A[0,N-1] = self.alpha_d1
        A[N-1,0] = self.alpha_d1
    
        B[0,N-1] = -1.
        B[N-1,0] = 1.
    
        
        Anew = A; Bnew = B;
        
    
    
        return Anew, Bnew


    def generic_difference(self, f, matrix, ngc):
        """
        ret = np.zeros(f.shape)
        ret[ngc:-ngc] = np.dot(matrix, f[ngc:-ngc])
        return ret
        """
        ret = np.dot(matrix, f)
        ret[0:ngc] = 0.
        ret[-ngc:] = 0.
        return ret


    def d_dx(self, f):
        """
        f is the function to differentiate 
        """
        if (self.grid.is_1d()):
            
            
            ngc = self.grid.ngc
            derivative = self.generic_difference(f, self.dx_matrix, ngc)

            derivative/=self.grid.dx
            
            if (self.grid.is_streched()):
                derivative = self.grid.d_eta_dx(derivative)
            
            
            return derivative
            
                
        if (self.grid.is_2d()):

            nx = self.grid.nx 
            ny = self.grid.ny
            ngc = self.grid.ngc
            derivative = np.zeros(f.shape)
            
            for j in range(ngc, ngc + ny):
                derivative[:, j] = self.generic_difference(f[:, j], self.dx_matrix, ngc)


            derivative/=self.grid.dx

            if (self.grid.stretched_x):
                derivative = self.grid.d_eta_dx(derivative)
            
            
            return derivative


    def d_dy(self, f):

        if (self.grid.is_1d()):
            print "Trying to compute d_dy on a 1D grid. Dont!"
            sys.exit()

        nx = self.grid.nx 
        ny = self.grid.ny
        ngc = self.grid.ngc
        derivative = np.zeros(f.shape)
        
        for i in range(ngc, ngc + nx):
            derivative[i,:] = self.generic_difference(f[i, :], self.dy_matrix, ngc)

        

        derivative/=(self.grid.dy)

        
        if (self.grid.is_streched()):
            pass
            #derivative = grid.d_eta_dy(derivative)
     
        
        return derivative

             
    def d_dxdx(self, f):
        """
        f is the function to differentiate 
        """
        if (self.grid.is_1d()):
            
            nx = self.grid.nx
            ngc = self.grid.ngc
            derivative = self.generic_difference(f, self.ddx_matrix, ngc)

            derivative/=(self.grid.dx*self.grid.dx)
            
            if (self.grid.is_streched()):
                derivative = self.grid.d_eta_dx(derivative)
            
            
            return derivative
            
                
        if (self.grid.is_2d()):

            nx = self.grid.nx 
            ny = self.grid.ny
            ngc = self.grid.ngc
            derivative = np.zeros(f.shape)
            
            for j in range(ngc, ngc + ny):
                derivative[:, j] = self.generic_difference(f[:,j], self.ddx_matrix, ngc)

            derivative/=(self.grid.dx*self.grid.dx)

            if (self.grid.is_streched()):
                pass
                #derivative = grid.d_eta_dx(derivative)
            
            
            return derivative        
            
        
            
        
    def d_dydy(self, f):

        if (self.grid.is_1d()):
            print "Trying to compute d_dydy on a 1D grid. Dont!"
            sys.exit()

        nx = self.grid.nx 
        ny = self.grid.ny
        ngc = self.grid.ngc
        derivative = np.zeros(f.shape)
        
        for i in range(ngc, ngc + nx):
            derivative[i,:] = self.generic_difference(f[i,:], self.ddy_matrix, ngc)

        derivative/=(self.grid.dy*self.grid.dy)

        
        if (self.grid.is_streched()):
            pass
            #derivative = grid.d_eta_dy(derivative)
     
        
        return derivative
    
    
            

    
     
                        
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
                
