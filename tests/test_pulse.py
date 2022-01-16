import numpy as np
import pytest
import copy
from mrpypulse import pulse
import sys
sys.path.append('..')


def test_pulse_init():

    p1 = pulse.Random()

    # assert np/tres/tp produce same output
    tp = p1.ns * p1.tres
    p2 = pulse.Pulse(x=p1.x, y=p1.y, tp=tp, tres=p1.tres)
    p3 = pulse.Pulse(x=p1.x, y=p1.y, ns=p1.ns, tp=tp)
    assert p1 == p2 == p3

    with pytest.raises(TypeError):
        p4 = pulse.Pulse(x=p1.x, y=p1.y, tres=p1.tres, ns=p1.ns, tp=p1.tp)

    with pytest.raises(TypeError):
        p4 = pulse.Pulse(x=p1.x, y=p1.y, tres=p1.tres)

    # assert x/y and r/ph produce same output
    r = np.sqrt(p1.x**2 + p1.y**2)
    ph = np.arctan2(p1.y, p1.x)
    p4 = pulse.Pulse(r=r, ph=ph, ns=p1.ns, tres=p1.tres)
    assert p1 == p4

    with pytest.raises(TypeError):
        p5 = pulse.Pulse(x=p1.x, y=p1.y, r=p1.r, ns=p1.ns, tres=p1.tres)

    with pytest.raises(TypeError):
        p5 = pulse.Pulse(x=p1.x, ns=p1.ns, tres=p1.tres)


def test_pulse_setattr():
    p1 = pulse.Random()
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


def test_pulse_add():
    p1 = pulse.Random()
    p2 = copy.deepcopy(p1)
    p3 = p1 + p2
    assert np.allclose(p3.x, p1.x+p2.x, rtol=1e-6, atol=1e-15)
    assert not np.isclose(p3.w1, p1.w1, rtol=1e-6, atol=1e-15)
    p1 += p2
    assert p1 == p3


def test_pulse_sub():
    p1 = pulse.Random()
    p2 = copy.deepcopy(p1)
    p3 = p1 - p2
    assert np.allclose(p3.x, p1.x-p2.x, rtol=1e-6, atol=1e-15)
    assert not np.isclose(p3.w1, p1.w1, rtol=1e-6, atol=1e-15)
    p1 -= p2
    assert p1 == p3


def test_pulse_add_ph_polyfit():
    p1 = pulse.Random()
    p2 = copy.deepcopy(p1)
    scale_ph = 2*np.random.rand()
    ph = np.random.uniform(low=-1.0, high=1.0, size=int(scale_ph * p1.ns))
    ph_corr = p1.add_ph_polyfit(ph)
    p2.ph += ph_corr
    assert p1 == p2


def test_hard_init():
    p1 = pulse.Hard(14e-6, 30e3)
    p2 = pulse.Pulse(tp=14e-6, ns=2,
                     r=np.array([30e3, 30e3]), ph=np.array([0, 0]))
    assert p1 == p2

    p1.w1 = 15e3
    assert p1 != p2


def test_parametrized_init():
    # assert chirped pulses Q, w1, tp, bw equivalence
    (am, fm, tr, tp, Q, bw) = ("WURST", "chirp", 1e-6, 500e-6, 5, 300e3)
    p1 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tr)
    p2 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, w1=p1.w1, tres=tr)
    p3 = pulse.Parametrized(AM=am, FM=fm, tp=tp, bw=bw, w1=p1.w1, tres=tr)
    p4 = pulse.Parametrized(AM=am, FM=fm, Q=Q, bw=bw, w1=p1.w1, tres=tr)
    assert p1 == p2 == p3 == p4

    # assert AM only and FM only pulses
    p5 = pulse.Parametrized(AM=None, FM=fm, Q=Q, bw=bw, w1=p1.w1, tres=tr)
    p6 = pulse.Parametrized(AM=am, FM=None, tp=tp, w1=p1.w1, tres=tr)
    p5.r = p6.r
    p6.ph = p5.ph
    assert p5 == p6

    # assert HS pulses Q, w1, tp, bw equivalence
    (am, fm, B) = ("sech", "sech", 10.6)
    p7 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tr, B=B)
    p8 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, w1=p7.w1, tres=tr, B=B)
    p9 = pulse.Parametrized(AM=am, FM=fm, tp=tp, bw=bw, w1=p7.w1, tres=tr, B=B)
    p10 = pulse.Parametrized(AM=am, FM=fm, Q=Q, bw=bw, w1=p7.w1, tres=tr, B=B)
    assert p7 == p8 == p9 == p10

    # smoothing percentage estimation
    (am, fm, tr) = ("superGaussian", "chirp", 0.25e-6)
    p11 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tr)
    assert p11.sm == 10.3
    (am, fm) = ("WURST", "chirp")
    p12 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tr)
    assert p12.sm == 10.7

    with pytest.raises(TypeError):
        p10 = pulse.Parametrized(AM=None, FM=None, tp=tp, w1=p1.w1, tres=tr)


def test_parametrized_setattr():
    # assert delta_t change
    (am, fm, tr, tp, Q, bw) = \
        ("Gaussian", "chirp", 2e-6, 480e-6, 3.7, 356.5e3)
    p1 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tr)
    p2 = copy.deepcopy(p1)
    p3 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tr,
                            start=p1.start + 500e-6)
    p1.delta_t += 500e-6
    assert p1.start == 500e-6 and p1 == p3
    p1.delta_t -= 250e-6
    p1.delta_t -= 250e-6
    assert p1 == p2

    # assert Q change
    p1.Q += 1.3
    p3 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=p1.Q, bw=bw, tres=tr)
    assert p1 == p3 and p1.w1 == p3.w1
    p1.Q -= 1.3
    assert p1 == p2 and p1.w1 == p2.w1

    # assert bw change
    p1.bw += 50e3
    p3 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=p1.bw, tres=tr)
    assert p1 == p3 and p1.w1 == p3.w1
    p1.bw -= 50e3
    assert p1 == p2 and p1.w1 == p2.w1

    # assert w1 change
    p1.w1 += 5e3
    p3 = pulse.Parametrized(AM=am, FM=fm, tp=tp, w1=p1.w1, bw=bw, tres=tr)
    assert p1 == p3 and p1.Q == p3.Q
    p1.w1 -= 5e3
    assert p1 == p2 and p1.Q == p2.Q

    # assert p change
    p1.p += 1
    p3 = pulse.Parametrized(AM=am, FM=fm, tp=tp, w1=p1.w1, bw=bw, tres=tr,
                            p=p1.p)
    assert p1 == p3
    p1.p -= 1
    assert p1 == p2

    # assert n change
    am = "WURST"
    p1 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tr)
    p2 = copy.deepcopy(p1)
    p1.n += 16
    p3 = pulse.Parametrized(AM=am, FM=fm, tp=tp, w1=p1.w1, bw=bw, tres=tr,
                            n=p1.n)
    assert p1 == p3 and p1.sm == p3.sm
    p1.n -= 16
    assert p1 == p2 and p1.sm == p2.sm
    am = "superGaussian"
    p1 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tr)
    p2 = copy.deepcopy(p1)
    p1.n += 16
    p3 = pulse.Parametrized(AM=am, FM=fm, tp=tp, w1=p1.w1, bw=bw, tres=tr,
                            n=p1.n)
    assert p1 == p3 and p1.sm == p3.sm
    p1.n -= 16
    assert p1 == p2 and p1.sm == p2.sm

    # assert sm change
    am = "sinsmoothed"
    p1 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tr)
    p2 = copy.deepcopy(p1)
    p1.sm += 5
    p3 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tr,
                            sm=p1.sm)
    assert p1 == p3
    p1.sm -= 5
    assert p1 == p2

    # assert B change
    (am, fm) = ("sech", "sech")
    p1 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tr)
    p2 = copy.deepcopy(p1)
    p1.B += 2.4
    p3 = pulse.Parametrized(AM=am, FM=fm, tp=tp, w1=p1.w1, bw=bw, tres=tr,
                            B=p1.B)
    assert p1 == p3 and p1.Q == p3.Q
    p1.B -= 2.4
    assert p1 == p2  # and p1.Q == p2.Q  # close but fails due - w1 change TODO

    # assert delta_f change
    p1.delta_f += 500
    p3 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tr,
                            delta_f=500)
    assert p1 == p3
    p1.delta_f -= 250
    p1.delta_f -= 250
    assert p1 == p2


def test_parametrized_str():
    (am, fm, tres, tp, Q, bw) = \
        ("sinsmoothed", "chirp", 1e-6, 856e-6, 6.74, 857.3e3)
    p1 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tres)
    p2 = copy.deepcopy(p1)
    assert str(p1) == str(p2)


def test_parametrized_reverse_sweep():
    (am, fm, tres, tp, Q, bw) = \
        ("superGaussian", "chirp", 0.3e-6, 400e-6, 2.61, 1152e3)
    p1 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tres)
    p2 = copy.deepcopy(p1)
    p1.reverse_sweep()
    assert np.allclose(p1.y, -p2.y, rtol=1e-6, atol=1e-15)
    assert not np.allclose(p1.ph, -p2.ph, rtol=1e-6, atol=1e-15)
    p1.reverse_sweep()
    assert p1 == p2


def test_parametrized_add_ph_polyfit():
    (am, fm, tres, tp, Q, bw) = \
        ("sinsmoothed", "chirp", 0.2e-6, 417e-6, 0.512, 470e3)
    p1 = pulse.Parametrized(AM=am, FM=fm, tp=tp, Q=Q, bw=bw, tres=tres)
    p2 = copy.deepcopy(p1)

    ns_ph = np.random.randint(20, 1000)
    ph = np.random.uniform(low=-1.0, high=1.0, size=ns_ph)
    ph_corr = p1.add_ph_polyfit(ph)
    p2.ph += ph_corr
    assert p1 == p2
