
import numpy as np
import matplotlib.pyplot as plt
import pulse as pulse

# hard pulse
p = pulse.Pulse(phi0=np.pi/2, tp=14e-6, ns=2, r=np.array([15e3, 15e3]), ph=np.array([0, 0]),ID=4)

# HS pulses
p2 = pulse.Parametrized(AM="tanh", FM="sech", tp=500e-6, Q=5., bw=300e3, tres=0.5e-6)
p3 = pulse.Parametrized(AM="tanh", FM="sech", tp=500e-6, Q=5., bw=300e3, ns=1000)
p2.phi0 = np.pi/2

ns = np.random.randint(2,1000)
x = np.random.uniform(low=-1.0, high=1.0, size=ns)
y = np.random.uniform(low=-1.0, high=1.0, size=ns)
tres = np.random.rand()
#p4 = pulse.Pulse(x=x, y=y, ns=ns, tres=tres)

tp = ns * tres
p4 = pulse.Pulse(x=x, y=y, tp=tp, tres=tres)
p5 = pulse.Pulse(x=x, y=y, ns=ns, tp=tp)

p4.plot(form='polar', title='Polar coordinates')

plt.figure()
plt.subplot(1,2,1)
p4.plot()
plt.subplot(1,2,2)
p5.plot(form='Cartesian')
