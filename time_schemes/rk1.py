import numpy as np


class Rk1:

    def __init__(self):
        self.type = "RK-1"
        

    def d_dt(self, f, df, dt):

        k1 = df(f)
        return dt*k1
