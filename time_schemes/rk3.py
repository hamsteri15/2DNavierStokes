import numpy as np


class Rk3:

    def __init__(self):
        self.type = "RK-3"
        

    def d_dt(self, f, df, dt):

        k1 = df(f)
        k2 = df(f + 0.5*k1*dt)
        k3 = df(f + (3./4.)*k2*dt)
        

        return (dt/9.)*(2*k1 + 3*k2 + 4*k3) 
