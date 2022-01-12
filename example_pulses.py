import numpy as np
import matplotlib.pyplot as plt
from mrpypulse import pulse
import matlab.engine
import copy

if __name__ == "__main__":

    # hard pulse
    p1 = pulse.Pulse(tp=14e-6, ns=2,
                     r=np.array([15e3, 15e3]), ph=np.array([0, 0]),
                     ID="hard 90­°", phi0=np.pi/2, )

    # random pulse
    ns = np.random.randint(2, 1000)
    x = np.random.uniform(low=-1.0, high=1.0, size=ns)
    y = np.random.uniform(low=-1.0, high=1.0, size=ns)
    tres = np.random.rand()
    tp = ns * tres
    p4 = pulse.Pulse(x=x, y=y, tp=tp, tres=tres, ID="random pulse")
    plt.figure()
    p4.plot(title=p4.ID)

    # chirp pulse
    p2 = pulse.Parametrized(AM="WURST", FM="chirp",
                            tp=100e-9, Q=4., bw=500e6, tres=0.25e-9)
    p2.phi0 = np.pi/2
    plt.figure()
    plt.subplot(1, 2, 1)
    p2.plot(form='polar')
    plt.subplot(1, 2, 2)
    p2.plot(form='Cartesian')
    plt.suptitle('WURST pulse')

    # HS pulses
    p3 = pulse.Parametrized(AM="sech", FM="sech",
                            tp=500e-6, Q=5., bw=300e3, ns=1000)
    plt.figure()
    p3.plot(title="HS pulse")

    # ESR Gaussian pulse
    p4 = pulse.Parametrized(AM="Gaussian", FM=None,
                            w1=40e6, tp=90e-9, tres=0.625e-9, p=15)
    plt.figure()
    p4.plot(form='polar')

    # ESR resonator compensated pulse (EasySpin)
    res_prof = np.genfromtxt('res_profile.txt',
                             delimiter='    ', skip_header=1)
    Hf = res_prof[:, 0]
    f = res_prof[:, 1]

    plt.figure()
    plt.plot(f, Hf)

    eng = matlab.engine.start_matlab()
    p2.resonator_easyspin(eng, f, Hf, 9.55)
    eng.quit()
    plt.figure()
    p2.plot()

    p1 = pulse.Random()
    p2 = copy.deepcopy(p1)
    p3 = p1 - p2

    plt.figure()
    p1.plot()
    plt.figure()
    p2.plot()
    plt.figure()
    p3.plot()
