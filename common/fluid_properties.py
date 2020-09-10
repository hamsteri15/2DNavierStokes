


class FluidProperties:

    def __init__(self, gamma, mu0, T0, Pr, R):
    
        self.gamma = gamma
        self.mu0 = mu0
        self.T0 = T0
        self.Pr = Pr
        self.R = R 
        
class FluidPropertiesAir(FluidProperties):

    def __init__(self):
        
        FluidProperties.__init__(self, 1.4, 273.1, 110.5, 0.71, 8.314)
        
