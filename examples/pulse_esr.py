import numpy as np
import matplotlib.pyplot as plt
from mrpypulse import pulse
from mrpypulse.magnetization import simulate

# ESR Gaussian pulse
p1 = pulse.Parametrized(AM="Gaussian", FM=None, w1=12e6,
                        tp=90e-9, tres=0.625e-9, p=15)

plt.figure()
p1.plot()

# ESR HS1 pulse
delta_f = -65e6
p2 = pulse.Parametrized(AM="tanh", FM="sech", w1=50e6, delta_f=delta_f,
                        bw=60e6, tp=100e-9, tres=0.625e-9, B=6/100e-9)

plt.figure()
p2.plot()

# simulation and plot of simulation results
off = np.linspace(delta_f-100e6, 50e6, 200)
magnHS = simulate([p2], off=off)
magnG = simulate([p1], off=off)

plt.figure()
plt.plot(off, magnHS[2, :], off, magnG[2, :])
