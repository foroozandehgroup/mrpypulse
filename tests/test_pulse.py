import numpy as np
import pytest
import copy
from mrpypulse import pulse
import sys
sys.path.append('..')


def random_pulse():
    """Generates a random pulse
    Returns
    -------
    p: pulse object
        random pulse
    """
    ns = np.random.randint(2, 1000)
    x = np.random.uniform(low=-1.0, high=1.0, size=ns)
    y = np.random.uniform(low=-1.0, high=1.0, size=ns)
    tres = np.random.rand()
    return pulse.Pulse(x=x, y=y, ns=ns, tres=tres)


def test_pulse_init():
    # random pulse parameters
    ns = np.random.randint(2, 1000)
    x = np.random.uniform(low=-1.0, high=1.0, size=ns)
    y = np.random.uniform(low=-1.0, high=1.0, size=ns)
    tres = np.random.rand()

    p1 = pulse.Pulse(x=x, y=y, ns=ns, tres=tres)

    # assert np/tres/tp produce same output
    tp = ns * tres
    p2 = pulse.Pulse(x=x, y=y, tp=tp, tres=tres)
    p3 = pulse.Pulse(x=x, y=y, ns=ns, tp=tp)
    assert p1 == p2 == p3

    with pytest.raises(TypeError):
        p4 = pulse.Pulse(x=x, y=y, tres=tres, ns=ns, tp=tp)

    with pytest.raises(TypeError):
        p4 = pulse.Pulse(x=x, y=y, tres=tres)

    # assert x/y and r/ph produce same output
    r = np.sqrt(x**2 + y**2)
    ph = np.arctan2(y, x)
    p4 = pulse.Pulse(r=r, ph=ph, ns=ns, tres=tres)
    assert p1 == p4

    with pytest.raises(TypeError):
        p5 = pulse.Pulse(x=x, y=y, r=r, ns=ns, tres=tres)

    with pytest.raises(TypeError):
        p5 = pulse.Pulse(x=x, ns=ns, tres=tres)


def test_pulse_setattr():
    p1 = random_pulse()
    p2 = copy.deepcopy(p1)

    # assert phi0 change
    p2.phi0 = np.random.uniform(low=-np.pi, high=np.pi)
    p1.ph += p2.phi0
    assert p1 == p2

    # assert w1 change
    p2.w1 = np.random.uniform(low=0.1, high=3)
    p1.r = p2.w1 * p1.r / np.max(p1.r)
    assert p1 == p2

    # assert start, end and t change
    position_shift = np.random.randint(1, 1000) * p1.tres
    p1.start = p1.start + position_shift
    p2.end = p2.end + position_shift
    assert np.allclose(p1.t, p2.t, rtol=1e-6, atol=1e-12)
    assert np.isclose(p1.start, p2.start, rtol=1e-6, atol=1e-12)
    assert np.isclose(p1.end, p2.end, rtol=1e-6, atol=1e-12)
    p1.t -= position_shift
    p2.start -= position_shift
    assert np.allclose(p1.t, p2.t, rtol=1e-6, atol=1e-12)
    assert np.isclose(p1.start, p2.start, rtol=1e-6, atol=1e-12)
    assert np.isclose(p1.end, p2.end, rtol=1e-6, atol=1e-12)

    # assert change of coordinates (x, y, r, ph)
    p1.x = p1.x * np.random.uniform(low=-1, high=1)
    p1.y = p1.y * np.random.uniform(low=-1, high=1)
    assert p1 != p2
    p1.r = p2.r
    p1.ph = p2.ph
    assert p1 == p2


def test_pulse_radd():
    p1 = random_pulse()
    p2 = copy.deepcopy(p1)
    p3 = p1 + p2
    assert np.allclose(p3.x, p1.x+p2.x, rtol=1e-6, atol=1e-15)
    assert not np.isclose(p3.w1, p1.w1, rtol=1e-6, atol=1e-15)
    p1 += p2
    assert p1 == p3


def test_pulse_add_ph_polyfit():
    p1 = random_pulse()
    p2 = copy.deepcopy(p1)
    scale_ph = 2*np.random.rand()
    ph = np.random.uniform(low=-1.0, high=1.0, size=int(scale_ph * p1.ns))
    ph_corr = p1.add_ph_polyfit(ph)
    p2.ph += ph_corr
    assert p1 == p2


def test_parametrized_init():
    # assert Q, w1, tp, bw equivalence
    (am, fm, tres, tp, Q, bw) = ("WURST", "chirp", 0.5e-6, 500e-6, 5, 300e3)
    p1 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tres)
    p2 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, w1=p1.w1, tres=tres)
    p3 = pulse.Parametrized(AM=am, FM=fm, tp=tp, bw=bw, w1=p1.w1, tres=tres)
    p4 = pulse.Parametrized(AM=am, FM=fm, Q=Q, bw=bw, w1=p1.w1, tres=tres)
    assert p1 == p2 == p3 == p4

    # assert AM only and FM only pulses
    p5 = pulse.Parametrized(AM=None, FM=fm, Q=Q, bw=bw, w1=p1.w1, tres=tres)
    p6 = pulse.Parametrized(AM=am, FM=None, tp=tp, w1=p1.w1, tres=tres)
    p5.r = p6.r
    p6.ph = p5.ph
    assert p5 == p6

    with pytest.raises(TypeError):
        p6 = pulse.Parametrized(AM=None, FM=None, tp=tp, w1=p1.w1, tres=tres)


def test_parametrized_str():
    (am, fm, tres, tp, Q, bw) = \
        ("sinsmoothed", "chirp", 0.5e-6, 500e-6, 5, 300e3)
    p1 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tres)
    p2 = copy.deepcopy(p1)
    assert str(p1) == str(p2)


def test_parametrized_reverse_sweep():
    (am, fm, tres, tp, Q, bw) = \
        ("superGaussian", "chirp", 0.5e-6, 500e-6, 5, 300e3)
    p1 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tres)
    p2 = copy.deepcopy(p1)
    p1.reverse_sweep()
    assert np.allclose(p1.y, -p2.y, rtol=1e-6, atol=1e-15)
    assert not np.allclose(p1.ph, -p2.ph, rtol=1e-6, atol=1e-15)
    p1.reverse_sweep()
    assert p1 == p2


def test_parametrized_add_ph_polyfit():
    (am, fm, tres, tp, Q, bw) = \
        ("sinsmoothed", "chirp", 0.5e-6, 500e-6, 5, 300e3)
    p1 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tres)
    p2 = copy.deepcopy(p1)

    ns_ph = np.random.randint(20, 1000)
    ph = np.random.uniform(low=-1.0, high=1.0, size=ns_ph)
    ph_corr = p1.add_ph_polyfit(ph)
    p2.ph += ph_corr
    assert p1 == p2
