""" SemiLagrangian.py : Implements semi-lagrangian integration schemes """
__author__ = "christian512"


import numpy as np

class SemiLagrangian1D:
    """
    SemiLagrangian scheme for periodic 1D problems.
    """
    def __init__(self,x,y,u,dt):
        """
        Initizalizes the SemiLagrangian object given initial configuration
        x : x positions (must be equidistant)
        y : variable at all x positions at time t
        u : speed of advection
        dt: time step
        """
        # Check dimensions
        assert x.shape[0] == y.shape[0], "x and y have different length"
        # Get the x stepsize
        dx = x[1] - x[0]
        # check if x is equidistant
        assert x[-1] == x[0] + dx*(x.shape[0]-1), "x array not equidistant"
        # Set class variables
        self.x = x
        self.y = y
        self.u = u
        self.dt = dt
        self.dx = dx

    def linear_evolve(self,nt=1):
        """
        Propagate the current variable using a linear interpolation
        nt : Number of time steps to be done
        """
        # loop through time steps
        for l in range(nt):
            # temporary array
            y_temp = np.empty(self.y.shape[0])
            # loop through variable array
            for i in range(self.y.shape[0]):
                # idx left to the departure point
                j = int(np.floor((self.x[i]-self.u*self.dt)/self.dx))
                # idx right to the departure point
                k = j+1
                # linear interpolation
                alpha = (self.x[i]-self.u*self.dt - j*self.dx)/self.dx
                y_temp[i] = (1-alpha)*self.y[j] + alpha*self.y[k]
            # copy array to current time
            self.y = np.copy(y_temp)
        #return current varibale
        return self.y

    def cubic_evolve(self,nt=1):
        """
        Propagates the variable using a cubic interpolation scheme.
        nt: number of time stepsize
        """
        #loop through time steps
        for l in range(nt):
            # temporary array
            y_temp = np.zeros(self.y.shape[0])
            # loop through array
            for i in range(self.y.shape[0]):
                # idx left to departure point
                x_dep = self.x[i]-self.u*self.dt
                j = int(np.floor(x_dep/self.dx))
                # alpha
                a = (self.x[i]-self.u*self.dt - j*self.dx)/self.dx
                # calculate next time step
                f = lambda x: x % self.y.shape[0] if x >= self.y.shape[0] else x
                y_temp[i] = - a * (1-a)*(2-a)/6 * self.y[f(j-1)]
                y_temp[i] += (1-a**2)*(2-a)/2 * self.y[f(j)]
                y_temp[i] += a*(1+a)*(2-a)/2 * self.y[f(j+1)]
                y_temp[i] -= a*(1-a**2)/6 * self.y[f(j+2)]
            self.y = np.copy(y_temp)
        return self.y
