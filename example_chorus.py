# script to generate the standard published CHORUS for NMR (2019)

import sequence as seq
import time

t90min = 500e-6
t180min = 1000e-6
bw = 300e3
tres = 0.5e-6
start = time.time()
chorus = seq.Exc_3fs(t90min, t180min, bw, tres, 
                     plot=True, t_del=0, polyfit=True)
print(f"Computation time: {format(time.time() - start,'0.2f')}s")
