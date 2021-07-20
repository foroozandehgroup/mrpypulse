import numpy as np
import pytest
#from mrpypulse.pulse import Pulse
import sys
sys.path.append('..')
import pulse

def test_pulse_init():
    """
    """
    # random pulse
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

    with pytest.raises(ValueError):
        p4=pulse.Pulse(x=x, y=y, tres=tres, ns=ns, tp=tp)

    with pytest.raises(ValueError):
        p4=pulse.Pulse(x=x, y=y, tres=tres)

    # assert x/y and r/ph produce same output
    r = np.sqrt(x**2 + y**2)
    ph = np.arctan2(y, x)
    p4 = pulse.Pulse(r=r, ph=ph, ns=ns, tres=tres)
    assert p1==p4

    with pytest.raises(ValueError):
        p5=pulse.Pulse(x=x, y=y, r=r, ns=ns, tres=tres)
    
    with pytest.raises(ValueError):
        p5=pulse.Pulse(x=x, ns=ns, tres=tres)
    

def pulse_setattr_test():
    """
    """
    # assert change in attribute values