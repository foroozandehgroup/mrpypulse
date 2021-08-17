import numpy as np
import matplotlib.pyplot as plt

import sys
  
# setting path
sys.path.append('../mrpypulse')

# importing
from mrpypulse.pulse import pulse

if __name__ == "__main__":

    # hard pulse
    p1 = pulse.Pulse(tp=14e-6, ns=2, r=np.array([15e3, 15e3]), ph=np.array([0, 0]),
                     ID="hard 90­°", phi0=np.pi/2, )
    
    # random pulse
    ns = np.random.randint(2,1000)
    x = np.random.uniform(low=-1.0, high=1.0, size=ns)
    y = np.random.uniform(low=-1.0, high=1.0, size=ns)
    tres = np.random.rand()
    tp = ns * tres
    p4 = pulse.Pulse(x=x, y=y, tp=tp, tres=tres, ID="random pulse")
    plt.figure()
    p4.plot(title=p4.ID)
    
    # chirp pulse
    p2 = pulse.Parametrized(AM="WURST", FM="chirp", 
                            tp=500e-6, Q=5., bw=300e3, tres=0.5e-6)
    p2.phi0 = np.pi/2
    plt.figure()
    plt.subplot(1,2,1)
    p2.plot(form='polar')
    plt.subplot(1,2,2)
    p2.plot(form='Cartesian')
    plt.suptitle('WURST pulse')
    
    # HS pulses
    p3 = pulse.Parametrized(AM="tanh", FM="sech", 
                            tp=500e-6, Q=5., bw=300e3, ns=1000)
    plt.figure()
    p3.plot(title="HS pulse")
