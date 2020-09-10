import numpy as np


class Rk4:

    def __init__(self):
        self.type = "RK-4"
        

    def d_dt(self, f, df):

        k1 = df(f)
        k2 = df(f + 0.5*k1*dt)
        k3 = df(f + 0.5*k2*dt)
        k4 = df(f + k3*dt)

        return (self.dt/6.)*(k1 + 2.*k2 + 2.*k3 + k4) 


    
