import unittest
import numpy as np
import sys
sys.path.append('../differentiators')
sys.path.append('../common')
from centered import CDN
from pade import Pade


class Test(unittest.TestCase):


    
    def test_1d(self):
        """
        runs 1d tests
        """

        print "---------1D TESTS----------"
        """
        print "Running scalar_1d_d_dx()"
        self.scalar_1d_d_dx()
        print "pass"

        print "Running vector_1d_d_dx()"
        self.vector_1d_d_dx()
        print "pass"

        print "Running scalar_1d_d_dxdx()"
        self.scalar_1d_d_dxdx()
        print "pass"

        print "Running vector_1d_d_dxdx()"
        self.vector_1d_d_dxdx()
        print "pass"
        """
    def test_2d(self):
        """
        runs 2d tests
        """
        print "---------2D TESTS----------"

        
        print "Running scalar_2d_d_dx()"
        self.scalar_2d_d_dx()
        print "pass"

        print "Running scalar_2d_d_dy()"
        self.scalar_2d_d_dy()
        print "pass"

        print "Running scalar_2d_d_dxdx()"
        self.scalar_2d_d_dxdx()
        print "pass"


        print "Running scalar_2d_d_dydy()"
        self.scalar_2d_d_dydy()
        print "pass"
        
        print "Running vector_2d_d_dx()"
        self.vector_2d_d_dx()
        print "pass"


        print "Running vector_2d_d_dy()"
        self.vector_2d_d_dy()
        print "pass"

        print "Running vector_2d_d_dxdx()"
        self.vector_2d_d_dxdx()
        print "pass"


        print "Running vector_2d_d_dydy()"
        self.vector_2d_d_dydy()
        print "pass"


    def set_up_1d(self):
        
        from grid_1d import Grid
        L = 10.0 #domain length good to make so that dx = 1
        nx = 10 #num points
        ngc = 5 #num ghost cells
        NX = nx + 2*ngc

        #differentiator objects to test
        #self.diffs = [CDN(6), Pade(nx)]
        self.grid = Grid(L, nx, ngc)

        #self.diffs = [CDN(self.grid, 2), CDN(self.grid, 4), CDN(self.grid, 6), Pade(self.grid)]
        self.diffs = [CDN(self.grid, 2), CDN(self.grid, 4), CDN(self.grid, 6)]

        
    
    def set_up_2d(self):
        
        

        from grid_2d import Grid
        L = 10.0 #domain length good to make so that dx = 1
        nx = 10 #num points
        ny = 10
        ngc = 5 #num ghost cells
        NX = nx + 2*ngc
        NY = ny + 2*ngc

        self.grid = Grid(L, nx, L, ny, ngc)

        #differentiators of order N
        #self.diffs = [CDN(self.grid, 2), CDN(self.grid, 4), CDN(self.grid, 6), Pade(self.grid)]
        self.diffs = [CDN(self.grid, 2), CDN(self.grid, 4), CDN(self.grid, 6)]

        

    

    def oned_array_tests(self, f, derivative, should_be_value, ngc, diff):

        #check that the return vector length is the same as input vector length
        self.assertTrue(derivative.size == f.size, msg="Scheme {0} returns a wrong length vector".format(diff.name))
        #check that ghost nodes are zero on both sides
        self.assertTrue(derivative[0:ngc].all() == 0.0, msg="Scheme {0} does not set left ghost cells to zero".format(diff.name))
        self.assertTrue(derivative[-ngc:].all() == 0.0, msg="order {0} does not set right ghost cells to zero".format(diff.name))

        #print derivative

        #sys.exit()

        #check that the return value is correct
        mean_derivative = np.mean(derivative[ngc:-ngc])
        msg = "{0} Does not differentiate correctly. Gives a mean of {1} it should be {2}".format(diff.name, mean_derivative, should_be_value)
        self.assertTrue(abs(mean_derivative - should_be_value) < 0.1, msg = msg)
        
    def twod_array_tests(self, f, derivative, should_be_value, ngc, diff):

        #check that the return vector length is the same as input vector length
        self.assertTrue(derivative.size == f.size, msg="Scheme {0} returns a wrong length vector".format(diff.name))
        #check that ghost nodes are zero on both sides
        #self.assertTrue(derivative[0:ngc].all() == 0.0, msg="Scheme {0} does not set left ghost cells to zero".format(diff.name))
        #self.assertTrue(derivative[-ngc:].all() == 0.0, msg="order {0} does not set right ghost cells to zero".format(diff.name))

        #check that the return value is correct
        mean_derivative = np.mean(derivative[ngc:-ngc, ngc:-ngc])
        msg = "{0} Does not differentiate correctly. Gives a mean of {1} it should be {2}".format(diff.name, mean_derivative, should_be_value)
        self.assertTrue(abs(mean_derivative - should_be_value) < 0.1, msg = msg)
        

    #############################
    #1D TESTS
    #############################

    def scalar_1d_d_dx(self):
        """
        Tests d_dx operator for an input vector of type f = [nx]
        """
        self.set_up_1d()
        ngc = self.grid.ngc

        f = np.arange(self.grid.NX)
        
        should_be_value = 1.0 #this is what the derivative should approximately be

        for diff in self.diffs:
            derivative = diff.d_dx(f)
            self.oned_array_tests(f, derivative, should_be_value, ngc, diff)


    def vector_1d_d_dx(self):
        """
        Tests d_dx operator for an input vector of type f = [nx, nq] i.e. 1d-vector field
        """
        nq = 5
        self.set_up_1d()
        ngc = self.grid.ngc
        f = np.zeros((self.grid.NX, nq))
        
        for i in range(self.grid.NX):
            for q in range(nq):
                f[i, q] = i
        
        should_be_value = 1.0 #this is what the derivative should approximately be

        for diff in self.diffs:
            derivative = diff.d_dx(f)
            self.oned_array_tests(f, derivative, should_be_value, ngc, diff)
            
        


    def scalar_1d_d_dxdx(self):
        """
        Tests d_dxdx operator for an input vector of type f = [nx]
        """
        self.set_up_1d()
        ngc = self.grid.ngc

        f = np.zeros(self.grid.NX)
        for i in range(self.grid.NX):
            f[i] = i*i
        
        should_be_value = 2.0 #this is what the derivative should approximately be

        for diff in self.diffs:
            derivative = diff.d_dxdx(f)
            self.oned_array_tests(f, derivative, should_be_value, ngc, diff)
            
        
       
    def vector_1d_d_dxdx(self):
        """
        Tests d_dxdx operator for an input vector of type f = [nx, nq] i.e. 1d-vector field
        """
        nq = 5
        self.set_up_1d()
        ngc = self.grid.ngc
        f = np.zeros((self.grid.NX, nq)) 
        
        for i in range(self.grid.NX):
            for q in range(nq):
                f[i, q] = i*i
        
        should_be_value = 2.0 #this is what the derivative should approximately be

        for diff in self.diffs:
            derivative = diff.d_dxdx(f)
            self.oned_array_tests(f, derivative, should_be_value, ngc, diff)
            
        
            

        
        
        

    #############################
    #2D TESTS
    #############################

    def scalar_2d_d_dx(self):
        """
        Tests d_dx operator for an input vector of type f = [nx,ny]
        """
        self.set_up_2d()
        ngc = self.grid.ngc
        
        f = np.zeros((self.grid.NX, self.grid.NY))

        
        for i in range(self.grid.NX):
            for j in range(self.grid.NY):
                f[i,j] = i
        
        should_be_value = 1.0

        for diff in self.diffs:
            derivative = diff.d_dx(f)
            self.twod_array_tests(f, derivative, should_be_value, ngc, diff)
            
        


    def scalar_2d_d_dxdx(self):
        """
        Tests d_dxdx operator for an input vector of type f = [nx,ny]
        """
        self.set_up_2d()
        ngc = self.grid.ngc
        
        f = np.zeros((self.grid.NX, self.grid.NY))
        
        for i in range(self.grid.NX):
            for j in range(self.grid.NY):
                f[i,j] = i*i
        
        should_be_value = 2.0

        for diff in self.diffs:
            derivative = diff.d_dxdx(f)
            self.twod_array_tests(f, derivative, should_be_value, ngc, diff)
            
        


    def scalar_2d_d_dy(self):
        """
        Tests d_dy operator for an input vector of type f = [nx,ny]
        """
        self.set_up_2d()
        ngc = self.grid.ngc
        
        f = np.zeros((self.grid.NX, self.grid.NY))

        for i in range(self.grid.NX):
            for j in range(self.grid.NY):
                f[i,j] = j

        should_be_value = 1.0

        for diff in self.diffs:
            derivative = diff.d_dy(f)
            self.twod_array_tests(f, derivative, should_be_value, ngc, diff)


    def scalar_2d_d_dydy(self):
        """
        Tests d_dy operator for an input vector of type f = [nx,ny]
        """
        self.set_up_2d()
        ngc = self.grid.ngc
        
        f = np.zeros((self.grid.NX, self.grid.NY))

        for i in range(self.grid.NX):
            for j in range(self.grid.NY):
                f[i,j] = j*j

        should_be_value = 2.0

        for diff in self.diffs:
            derivative = diff.d_dydy(f)
            self.twod_array_tests(f, derivative, should_be_value, ngc, diff)        
        

    def vector_2d_d_dx(self):
        """
        Tests d_dx operator for an input vector of type f = [nx,ny,nq]
        """
        nq = 5

        self.set_up_2d()
        ngc = self.grid.ngc

        f = np.zeros((self.grid.NX, self.grid.NY, nq))
        
        for i in range(self.grid.NX):
            for j in range(self.grid.NY):
                for q in range(nq):
                    f[i,j,q] = i

        should_be_value = 1.0

        for diff in self.diffs:
            derivative = diff.d_dx(f)
            self.twod_array_tests(f, derivative, should_be_value, ngc, diff)


    def vector_2d_d_dy(self):
        """
        Tests d_dy operator for an input vector of type f = [nx,ny,nq]
        """
        nq = 5

        self.set_up_2d()
        ngc = self.grid.ngc

        f = np.zeros((self.grid.NX, self.grid.NY, nq))
        
        for i in range(self.grid.NX):
            for j in range(self.grid.NY):
                for q in range(nq):
                    f[i,j,q] = j

        should_be_value = 1.0

        for diff in self.diffs:
            derivative = diff.d_dy(f)
            self.twod_array_tests(f, derivative, should_be_value, ngc, diff)



        


    def vector_2d_d_dxdx(self):
        """
        Tests d_dxdx operator for an input vector of type f = [nx,ny,nq]
        """
        nq = 5

        self.set_up_2d()
        ngc = self.grid.ngc

        f = np.zeros((self.grid.NX, self.grid.NY, nq))

        for i in range(self.grid.NX):
            for j in range(self.grid.NY):
                for q in range(nq):
                    f[i,j,q] = i*i #this is a parabola and the derivative should be close to 2.0

        should_be_value = 2.0

        for diff in self.diffs:
            derivative = diff.d_dxdx(f)
            self.twod_array_tests(f, derivative, should_be_value, ngc, diff)




    def vector_2d_d_dydy(self):
        """
        Tests d_dydy operator for an input vector of type f = [nx,ny,nq]
        """
        nq = 5

        self.set_up_2d()
        ngc = self.grid.ngc

        f = np.zeros((self.grid.NX, self.grid.NY, nq))

        for i in range(self.grid.NX):
            for j in range(self.grid.NY):
                for q in range(nq):
                    f[i,j,q] = j*j #this is a parabola and the derivative should be close to 2.0

        should_be_value = 2.0

        for diff in self.diffs:
            derivative = diff.d_dydy(f)
            self.twod_array_tests(f, derivative, should_be_value, ngc, diff)


    



    


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()    