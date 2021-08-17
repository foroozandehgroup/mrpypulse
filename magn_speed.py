
import numpy as np
import matplotlib.pyplot as plt
from mrpypulse import pulse
from mrpypulse.magnetization import simulate, simulate2, plot_magn
from mrpypulse import sequence as seq
import time

if __name__ == "__main__":

    # comparison between simulate and simulate2
    
    t90min = 100e-6
    t180min = 100e-6
    bw = 500e3
    tres = 2e-6
    chorus = seq.Exc_3fs(t90min, t180min, bw, tres, plot=False, polyfit=True)
    
    # check that magnetization is recocused in centre part
    off = np.linspace(-0.5*bw, 0.5*bw, 100)
    start = time.time()
    magn, off = simulate(chorus.pulses, offsets=off, pc=chorus.pc)
    print(f"simulate computation time: {format(time.time() - start,'0.2f')}s")
    plot_magn(magn, off)
    
    start = time.time()
    start = time.time()
    magn, off = simulate2(chorus.pulses, offsets=off, pc=chorus.pc)
    print(f"simulate computation time: {format(time.time() - start,'0.2f')}s")
    plot_magn(magn, off)
