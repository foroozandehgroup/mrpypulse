import numpy as np
import matplotlib.pyplot as plt
from mrpypulse import pulse
import matlab.engine

# ESR resonator compensated pulse (requires matlab.engine, MATLAB, EasySpin)
res_prof = np.genfromtxt('res_profile.txt',
                         delimiter='    ', skip_header=1)
Hf = res_prof[:, 0]
f = res_prof[:, 1]

plt.figure()
plt.plot(f, Hf)
plt.title('resonator profile')

# chirp pulse
p = pulse.Parametrized(AM="WURST", FM="chirp",
                       tp=100e-9, Q=4., bw=500e6, tres=0.25e-9)


eng = matlab.engine.start_matlab()
p.resonator_easyspin(eng, f, Hf, 9.55)
eng.quit()
plt.figure()
p.plot()
