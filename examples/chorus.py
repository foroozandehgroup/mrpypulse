# Script to generate the standard published CHORUS for NMR
# Foroozandeh, M., et al. (2019). "Improved ultra-broadband chirp excitation."
# Journal of Magnetic Resonance 302: 28-33.
#
# can take up to a few minutes to run.


from mrpypulse import sequence as seq
import time

t90min = 500e-6
t180min = 1000e-6
bw = 300e3
tres = 0.5e-6
start = time.time()
chorus = seq.Exc_3fs(t90min, t180min, bw, tres,
                     plot=True, t_del=0, polyfit=True,
                     pulse_args={"sm": 10},
                     polyfit_args={"deg": 10})

print(f"Computation time: {format(time.time() - start,'0.2f')}s")
