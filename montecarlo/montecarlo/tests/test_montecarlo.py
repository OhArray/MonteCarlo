"""
Unit and regression test for the montecarlo package.
"""

# Import package, test suite, and other packages as needed
import montecarlo
import pytest
import sys
import random
import numpy as np
import copy as cp

def test_montecarlo_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "montecarlo" in sys.modules


def testa():
    assert 1 == 1

def test_average_values():
    N=10
    conf = montecarlo.BitString([0]*N)
    #conf.initialize(M=20)

    T = 2.0

    
    # now test the general ising hamiltonian
    Jval = 1.0
    mu = [.1 for i in range(N)]
    J = []
    for i in range(N):
        J.append([((i+1) % N, Jval), ((i-1) % N, Jval)])
    ham2 = montecarlo.IsingHamiltonian(J=J, mu=mu)
    E, M, HC, MS = ham2.compute_average_values(conf, T) 
 
    assert(np.isclose(E, -4.6378514858094695))
    assert(np.isclose(M, -0.1838233606011354 ))
    assert(np.isclose(HC, 1.9883833749653714 ))
    assert(np.isclose(MS, 1.8391722085614428))


def test_metropolis():
    random.seed(2)
    N=20
    conf = montecarlo.BitString([0]*N)
    T = 2

    J = []
    Jval = 1.0
    mu = [.1 for i in range(N)]
    for i in range(N):
        J.append([((i+1) % N, Jval), ((i-1) % N, Jval)])
    ham = montecarlo.IsingHamiltonian(J=J, mu=mu)

    conf = montecarlo.BitString([0]*N)
    E, M, EE, MM = montecarlo.metropolis_montecarlo(ham, conf, T=T, nsweep=8000, nburn=1000)
    HC = (EE[-1] - E[-1]*E[-1])/T/T
    MS = (MM[-1] - M[-1]*M[-1])/T
    print
    print("     E:  %12.8f" %(E[-1]))
    print("     M:  %12.8f" %(M[-1]))
    print("     HC: %12.8f" %(HC))
    print("     MS: %12.8f" %(MS))

    assert(np.isclose(-9.31, E[-1]))
     

def test_classes():
    random.seed(2)
    conf = montecarlo.BitString([0]*10)
    conf.set_magnetization(M=5)
    assert(all(conf.array() == [1, 1, 1, 0, 0, 0, 0, 1, 1, 0]))

    N = 10
    J = []
    Jval = -1.0
    mu = [-.001 for i in range(N)]
    for i in range(N):
        J.append([((i+1) % N, Jval), ((i-1) % N, Jval)])
    ham2 = montecarlo.IsingHamiltonian(J=J, mu=mu)

    e = ham2.energy(conf)
    print(" Energy = ", e)
    assert(np.isclose(e,-2))

    conf.flip_site(3)
    print(conf.array())
    print(" Energy = ", e)
    assert(np.isclose(ham2.energy(conf),-2.002))
    
    # now flip back
    conf.flip_site(3)
    print(conf.array())
    e = ham2.energy(conf)
    print(" Energy = ", e)
    assert(np.isclose(ham2.energy(conf),-2.00))


    conf_old = cp.deepcopy(conf)
   
    ham2.mu = np.array([-1.1 for i in range(N)])
    ham2.metropolis_sweep(conf, T=.9)
    print(conf_old, " --> ", conf)
    print("Energy: %12.8f --> %12.8f" %(e, ham2.energy(conf)))  
    assert(all(conf.array() == np.ones(10)))
    

    random.seed(2)
    conf.set_int(44)
    conf_old = cp.deepcopy(conf)

    
    Jval = 1.0
    mu = [-.1 for i in range(N)]
    J = []
    for i in range(N):
        J.append([((i+1) % N, Jval), ((i-1) % N, Jval)])
    ham2 = montecarlo.IsingHamiltonian(J=J, mu=mu)
    random.seed(2)
    conf.set_int(44)
    conf_old = cp.deepcopy(conf)
    ham2.metropolis_sweep(conf, T=.9)
    print(conf_old, " --> ", conf)
    print(conf_old)
    print(ham2.energy(conf_old))  
    print("Energy: %12.8f --> %12.8f" %(e, ham2.energy(conf)))  
    assert(all(conf.array() == [1,1,1,0,1,0,0,1,1,0]))
   
def test_delta_e():
    random.seed(2)
    N = 20
    conf = montecarlo.BitString([0]*N)
    conf.set_magnetization(M=10)

    Jval = -1.0

    mu = [-.001 for i in range(N)]
    J = []
    for i in range(N):
        J.append([((i+1) % N, Jval), ((i-1) % N, Jval)])


    ham = montecarlo.IsingHamiltonian(J=J, mu=mu)

    e1 = ham.energy(conf)
    print(conf)
    print(" Energy = ", e1)

    delta_e1 = ham.delta_e_for_flip(3, conf)
    delta_e2 = ham.delta_e_for_flip_slow(3, conf)
    delta_e3 = montecarlo.delta_e_for_flip_fast(3, conf.array(), 
                                         ham.nodes[3],
                                         ham.js[3],
                                         ham.mu)

    conf.flip_site(3)
    e2 = ham.energy(conf)
    print(conf)
    print(" Energy = ", e2)

    print(" delta E: %12.8f" %(e2-e1))
    print(" delta E: %12.8f" %(delta_e1))
    print(" delta E: %12.8f" %(delta_e2))
    print(" delta E: %12.8f" %(delta_e3))

    assert(np.isclose(e2-e1, delta_e1))
    assert(np.isclose(e2-e1, delta_e2))
    assert(np.isclose(e2-e1, delta_e3))
    
if __name__== "__main__":
    test_montecarlo_imported()
    test_classes()
    test_average_values()
    test_metropolis()
    test_delta_e()
