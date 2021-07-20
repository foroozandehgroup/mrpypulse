import numpy as np
import pytest
import copy
import sys
sys.path.append('..')
import pulse

def random_pulse():
    """Generates a random pulse
    Returns
    -------
    p: pulse object
        random pulse
    """
    ns = np.random.randint(2,1000)
    x = np.random.uniform(low=-1.0, high=1.0, size=ns)
    y = np.random.uniform(low=-1.0, high=1.0, size=ns)
    tres = np.random.rand()
    return pulse.Pulse(x=x, y=y, ns=ns, tres=tres)

def test_pulse_init():
    # random pulse parameters
    ns = np.random.randint(2,1000)
    x = np.random.uniform(low=-1.0, high=1.0, size=ns)
    y = np.random.uniform(low=-1.0, high=1.0, size=ns)
    tres = np.random.rand()

    p1 = pulse.Pulse(x=x, y=y, ns=ns, tres=tres)

    # assert np/tres/tp produce same output
    tp = ns * tres
    p2 = pulse.Pulse(x=x, y=y, tp=tp, tres=tres)
    p3 = pulse.Pulse(x=x, y=y, ns=ns, tp=tp)
    assert p1==p2==p3

    with pytest.raises(TypeError):
        p4=pulse.Pulse(x=x, y=y, tres=tres, ns=ns, tp=tp)

    with pytest.raises(TypeError):
        p4=pulse.Pulse(x=x, y=y, tres=tres)

    # assert x/y and r/ph produce same output
    r = np.sqrt(x**2 + y**2)
    ph = np.arctan2(y, x)
    p4 = pulse.Pulse(r=r, ph=ph, ns=ns, tres=tres)
    assert p1==p4

    with pytest.raises(TypeError):
        p5=pulse.Pulse(x=x, y=y, r=r, ns=ns, tres=tres)
    
    with pytest.raises(TypeError):
        p5=pulse.Pulse(x=x, ns=ns, tres=tres)
    

def test_pulse_setattr():
    p1 = random_pulse()
    p2 = copy.deepcopy(p1)

    # assert phi0 change
    p2.phi0 = np.random.uniform(low=-np.pi, high=np.pi)
    p1.ph += p2.phi0
    assert p1==p2

    # assert w1 change
    p2.w1 = np.random.uniform(low=0.1, high=3)
    p1.r = p2.w1 * p1.r / np.max(p1.r)
    assert p1==p2

    # assert change of coordinates
    p1.x = p1.x * np.random.uniform(low=-1, high=1)
    p1.y = p1.y * np.random.uniform(low=-1, high=1)
    assert p1!=p2
    p1.r = p2.r
    p1.ph = p2.ph
    assert p1==p2