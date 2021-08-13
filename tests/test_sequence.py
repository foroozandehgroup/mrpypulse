import numpy as np
import simulation
import sequence

def test_Exc_3fs_init():
    """
    Generates a short chorus with a low time resolution and test that the 
    magnetization is grossly refocused on y.
    """
    t90min = 100e-6
    t180min = 101e-6
    bw = 506e3
    tres = 2e-6
    chorus = sequence.Exc_3fs(t90min, t180min, bw, tres, plot=False, polyfit=True)
    
    # check that magnetization is recocused in centre part
    off = np.linspace(-0.25*bw, 0.25*bw, 10)
    magn, off = simulation.simulate(chorus.pulses, offsets=off, pc=chorus.pc)
    assert np.all(magn[1,:] > 0.9)