# --------------------------------------------------------------------------- #
# Demonstration file on use of mrpypulse

import time
import numpy as np
import matplotlib.pyplot as plt
from mrpypulse import pulse, sequence
from mrpypulse.magnetization import simulate, plot_magn

# --------------------------------------------------------------------------- #
# Create pulses with pulse module

# creation of an HS1 pulse for inversion in NMR
p1 = pulse.Parametrized(AM='sech', FM='sech', bw=65e3,
                        tp=1000e-6, Q=0.441, tres=1e-6,
                        ID='demo_pulse')

# look at the pulse
plt.figure()
p1.plot()

# --------------------------------------------------------------------------- #
# Manipulate pulses with their methods

# modify the pulse adiabaticity factor
p1.Q = 5

plt.figure()
p1.plot()

# --------------------------------------------------------------------------- #
# Simulate pulses with functions from the magnetization module

# simulate the effect of a single pulse
offsets = np.linspace(-p1.bw/2, p1.bw/2, 50)
magnetization = simulate([p1], off=offsets)

plot_magn(magnetization, offsets)

# --------------------------------------------------------------------------- #
# Use several pulses to make sequences
# double pulse for refocusing

# new position
p2 = pulse.Parametrized(AM='sech', FM='sech', bw=65e3,
                        tp=1000e-6, Q=5, tres=1e-6,
                        start=p1.end)

plt.figure()
p2.plot()
plt.xlim(0, p2.end)

# --------------------------------------------------------------------------- #
# Create sequences directly with the sequence module

# ABSTRUSE sequence creation
abstruse = sequence.Exc_3fs(t90min=p1.tp/2, t180min=p1.tp, bw=p1.bw,
                            tres=p1.tres,
                            pulse_args={'AM': 'sech', 'FM': 'sech'})

# plot the sequence
plt.figure()
abstruse.plot()

# --------------------------------------------------------------------------- #
# Use additional options from the sequence creation functions

# polynomial fitting phase correction
start = time.time()
abstruse = sequence.Exc_3fs(t90min=p1.tp/2, t180min=p1.tp, bw=p1.bw,
                            tres=p1.tres, ID='demo_seq',
                            pulse_args={'AM': 'sech', 'FM': 'sech'},
                            polyfit=True, plot=True)

print(f"Computation time: {format(time.time() - start,'0.2f')}s")

# --------------------------------------------------------------------------- #
# Manipulate sequences with their methods

# add a pulse at the end of the sequence
abstruse.append(p1)

plt.figure()
abstruse.plot()

# --------------------------------------------------------------------------- #
# Export a pulse shape to be used on the instrument

p1.topspin_file()

# --------------------------------------------------------------------------- #
# Export a whole sequence

p1.ID = 'p4'
abstruse.seq2topspin()
