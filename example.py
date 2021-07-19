
import numpy as np
import matplotlib.pyplot as plt
import pulse as pulse
#p = pulse.Pulse(phi0=np.pi/2, tp=14e-6, ns=2, r=np.array([15e3, 15e3]), ph=np.array([0, 0]))


p2 = pulse.Parametrized(AM="tanh", FM="sech", tp=500e-6, Q=5, bw=300e3, tres=0.5e-6)

print(p2.w1)
p2.plot()

p2.w1 = 2*p2.w1
print(p2.w1)
p2.plot()

plt.show()