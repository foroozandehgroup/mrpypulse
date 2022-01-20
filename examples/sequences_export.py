# example showing how to export pulse sequences to Xepr and TopSpin shapes

from mrpypulse import sequence


# chorus for ESR
chorus = sequence.Exc_3fs(t90min=80e-9, t180min=160e-9, bw=450e6,
                          tres=0.625e-9, plot=False, polyfit=False,
                          ID='chorus_epr')

chorus.seq2xepr()

# chorus for NMR
chorus = sequence.Exc_3fs(t90min=500e-6, t180min=1000e-6, bw=500e3,
                          tres=0.5e-6, plot=False, polyfit=False,
                          ID='chorus_nmr')

chorus.seq2topspin()
