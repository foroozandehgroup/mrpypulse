import numpy as np
import copy
from mrpypulse import pulse, sequence, magnetization


def test_pc_rec():
    pc1 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    pc2 = np.array([0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3])
    pc3 = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3])

    ph_rec = sequence.pc_rec([pc1, pc2, pc3], [-1, +2, -2])

    pc31 = np.array([0, 2, 0, 2, 2, 0, 2, 0, 0, 2, 0, 2, 2, 0, 2, 0])

    assert np.all(ph_rec == pc31)

    pc1 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2])
    pc2 = np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1])
    pc3 = np.array([0, 0, 1, 1, 2, 2, 3, 3, 0, 0, 1, 1, 2, 2, 3, 3])

    ph_rec = sequence.pc_rec([pc1, pc2, pc3], [-1, +2, -2])

    pc31 = np.array([0, 2, 2, 0, 0, 2, 2, 0, 2, 0, 0, 2, 2, 0, 0, 2])

    assert np.all(ph_rec == pc31)


def test_sequence_insert():
    t90min = 80e-6
    t180min = 101e-6
    bw = 506e3
    tres = 0.5e-6
    chorus = sequence.Exc_3fs(t90min, t180min, bw, tres,
                              plot=False, polyfit=False, t_del=80e-6)
    chorus.insert(0, 100e-6)
    p = pulse.Parametrized(AM="sech", FM="sech",
                           tp=500e-6, Q=5., bw=300e3, ns=1000)
    chorus.insert(1, p)

    chorus2 = sequence.Exc_3fs(t90min, t180min, bw, tres,
                               plot=False, polyfit=False, t_del=0e-6)

    chorus2.insert(0, 100e-6)
    chorus2.insert(2, 80e-6)
    chorus2.insert(3, 80e-6)
    chorus2.insert(1, p)

    assert chorus == chorus2


def test_sequence_append():
    chorus = sequence.Exc_3fs(t90min=79e-6, t180min=180e-6, bw=487.52e3,
                              tres=0.5e-6, plot=False, polyfit=False)
    p = pulse.Parametrized(AM="sech", FM="sech",
                           tp=500e-6, Q=5., bw=300e3, ns=1000)
    chorus2 = copy.deepcopy(chorus)
    chorus2.insert(3, 100e-6)
    chorus2.insert(3, p)
    chorus.append(100e-6)
    chorus.append(p)

    assert chorus == chorus2


def test_Exc_3fs_init():
    """
    Generates a short chorus with a low time resolution and test that the
    magnetization is grossly refocused on y.
    """
    t90min = 80e-6
    t180min = 101e-6
    bw = 506e3
    tres = 2e-6
    chorus = sequence.Exc_3fs(t90min, t180min, bw, tres,
                              plot=False, polyfit=True)

    # check that magnetization is recocused in centre part
    off = np.linspace(-0.25*bw, 0.25*bw, 10)
    magn = magnetization.simulate(
        chorus.pulses, off=off, pc=chorus.pc)
    assert np.all(magn[1, :] > 0.9)
