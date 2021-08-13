
import numpy as np
import matplotlib.pyplot as plt
import pulse as pulse
import simulation
import sequence as seq
import time


# comparison between simulate and simulate2

t90min = 100e-6
t180min = 100e-6
bw = 500e3
tres = 2e-6
chorus = seq.Exc_3fs(t90min, t180min, bw, tres, plot=False, polyfit=True)

# check that magnetization is recocused in centre part
off = np.linspace(-0.5*bw, 0.5*bw, 100)
start = time.time()
magn, off = simulation.simulate(chorus.pulses, offsets=off, pc=chorus.pc)
print(f"simulate computation time: {format(time.time() - start,'0.2f')}s")
simulation.plot_magn(magn, off)

start = time.time()
start = time.time()
magn, off = simulation.simulate2(chorus.pulses, offsets=off, pc=chorus.pc)
print(f"simulate computation time: {format(time.time() - start,'0.2f')}s")
simulation.plot_magn(magn, off)
