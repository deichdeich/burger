import numpy as np
import matplotlib.pyplot as plt

class berger_solve(object):
    def __init__(self,
                 dx = 0.1,
                 dt = 0.01,
                 mu = 0.01,
                 tmax = 1):
        """
        simulation parameters
        """
        self.dx = dx
        self.dt = dt
        self.tmax = tmax
        self.nsteps = int(tmax / dt)
        self.mu = mu
        
        """
        the thing that appears when you solve for the next timestep
        """
        self.factor = (self.dt) / (2 * self.dx)
        
        """
        initializing the grid
        """
        self.grid = self.inish_condish()
        
        self.Nx = len(self.grid) # get the grid spacing
        
        self.clock = 0

    def _set_dt(self, new_dt):
    
        self.dt = new_dt
        self.factor = (self.dt) / (2 * self.dx)
        grid = self.inish_condish()
        self.nsteps = int(self.tmax / new_dt)
        self.clock = 0
        return(new_dt)
    
    def _set_dx(self, new_dx):
        self.dx = new_dx
        self.factor = (self.dt) / (2 * self.dx)
        grid = self.inish_condish()
        self.Nx = len(self.grid)
        self.clock = 0
        return(new_dt)

    """
    inish_condish returns a 1xNx array which holds the value of the field
    """
    def inish_condish(self):
        # these lines set up the initial condition
        grid = np.arange(0, 1, self.dx) 
        term1 = (1/2) * np.sin(np.pi * grid)
        term2 = np.sin(2 * np.pi * grid)
        grid = term1 + term2
        
        # create dummy cells at the edges to impose dirichlet b.c.'s
        grid = self.impose_dbc(grid)
        
        return(term1 + term2)
    
    """
    impose_dbc takes a grid and returns a grid with two cells at either end which match
    the end points of the input grid.  This imposes the dirichlet boundary conditions,
    assuming the initial condition satisfies those.
    """
    def impose_dbc(self, partial_grid):
        grid = np.zeros(len(partial_grid) + 2)
        grid[1:-1] = partial_grid
        grid[0] = partial_grid[0]
        grid[-1] = partial_grid[-1]
        return(grid)
    
    def update(self, grid):
        new_grid = np.zeros_like(grid)
        
        # this enforces the boundary conditions; the edges won't be touched.
        for i in range(1, self.Nx - 1):
            new_grid[i] = grid[i] + self.factor *\
                          (\
                           (\
                            (self.mu / self.dx) *\
                            (grid[i+1] - (2 * grid[i]) + grid[i-1])\
                           )-\
                          grid[i] * (grid[i+1] - grid[i-1])\
                          )
        return(new_grid)
    
    def evolve(self):
        history = np.zeros((self.nsteps, self.Nx))
        for step in range(self.nsteps):
            history[step] = self.grid
            self.grid = self.update(self.grid)
            self.clock += self.dt
        
        return(history)
            
        