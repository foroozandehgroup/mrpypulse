import pulse
import magnetization
import numpy as np


def test_simulate():
    # inversion pulse
    p = pulse.Parametrized(AM="sinsmoothed", FM="chirp", sm=20,
                           tp=250e-6, Q=5., bw=600e3, tres=2e-6)
    # simulation
    off = np.linspace(-p.bw/2, p.bw/2, 50)
    magn, offsets = magnetization.simulate([p],offsets=off)
    Mz = magn[2,:]
    
    # check that magnetization is inverted in centre part
    assert np.all(Mz[20:30] < -0.95)
