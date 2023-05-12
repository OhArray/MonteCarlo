from numba import jit

import numpy as np
import random
import montecarlo

@jit(nopython=True)
def delta_e_for_flip_fast(i, config, nodes, J, mu):
    """
    Parameters
    ----------
    i        : int
        Index of site to flip
    config   : list
        list of 0's and 1's 
    nodes    : [[]]]
        for each node, list of nodes connected
    J        : [[]]
        for each node, list of J couplings to connected nodes
    mu       : list
        for each node, strength of local field, mu
    
    Returns
    -------
    del_e  : float
        Returns the energy change
    """
    
    del_e = 0.0
    del_si = 2
    if config[i] == 1:
        del_si = -2

    for jidx, j in enumerate(nodes):
        del_e += (2.0*config[j]-1.0) * J[jidx] * del_si

    del_e += mu[i] * del_si 
    return del_e

class IsingHamiltonian:

    def __init__(self, J=[[()]], mu=np.zeros(1)):
        """ Constructor 
    
        Parameters
        ----------
        J: list of lists of tuples, optional
            Strength of coupling, e.g, 
            [(4, -1.1), (6, -.1)]
            [(5, -1.1), (7, -.1)]
        mu: vector, optional
            local fields 
        """
        self.J = J
        self.mu = mu
    
        self.nodes = []
        self.js = []
 
        for i in range(len(self.J)):
            self.nodes.append(np.zeros(len(self.J[i]), dtype=int))
            self.js.append(np.zeros(len(self.J[i])))
            for jidx,j in enumerate(self.J[i]):
                self.nodes[i][jidx] = j[0]
                self.js[i][jidx] = j[1]
        self.mu = np.array([i for i in self.mu])

    def energy(self, config):
        """Compute energy of configuration, `config` 
            
            .. math::
                E = \\left<\\hat{H}\\right> 

        Parameters
        ----------
        config   : BitString
            input configuration 

        Returns
        -------
        energy  : float
            Energy of the input configuration
        """
            
        e = 0.0

        for i in range(config.N):
            for j in self.J[i]:
                if j[0] < i:
                    continue
                if config[i] == config[j[0]]:
                    e += j[1]
                else:
                    e -= j[1]

        e += np.dot(self.mu, 2*config.array()-1)

        return e
    
    def compute_average_values(self, conf, T):
        """ Compute Average values exactly

        Parameters
        ----------
        conf   : :class:`BitString`
            input configuration 
        T      : int
            Temperature
        
        Returns
        -------
        E  : float 
            Energy
        M  : float
            Magnetization
        HC : float
            Heat Capacity
        MS : float
            Magnetic Susceptability
        """
        E  = 0.0
        M  = 0.0
        Z  = 0.0
        EE = 0.0
        MM = 0.0

        for i in range(conf.n_dim):
            conf.set_int(i)
            Ei = self.energy(conf)
            Zi = np.exp(-Ei/T)
            E += Ei*Zi
            EE += Ei*Ei*Zi
            Mi = np.sum(2*conf.array()-1)
            M += Mi*Zi
            MM += Mi*Mi*Zi
            Z += Zi
        
        E = E/Z
        M = M/Z
        EE = EE/Z
        MM = MM/Z
        
        HC = (EE - E*E)/(T*T)
        MS = (MM - M*M)/T
        return E, M, HC, MS
    
    
    def metropolis_sweep(self, conf, T=1.0):
        """Perform a single sweep through all the sites and return updated configuration

        Parameters
        ----------
        conf   : :class:`BitString`
            input configuration 
        T      : int
            Temperature
        
        Returns
        -------
        conf  : :class:`BitString`
            Returns updated config
        """
 
        for site_i in range(conf.N):

            delta_e = delta_e_for_flip_fast(site_i, conf.array(), 
                                            self.nodes[site_i], 
                                            self.js[site_i], 
                                            self.mu)      

            accept = True
            if delta_e > 0.0:
                rand_comp = random.random()
                if rand_comp > np.exp(-delta_e/T):
                    accept = False
            if accept:
                conf.flip(site_i)
        return conf