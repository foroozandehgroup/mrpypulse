
"This module contains functions to simulate pulses and sequences."

import numpy as np

# TODO code this module

def simulate(self, off, time=None, B1=None, use_pc=False, magn_init=None):
    """Simulation of the magnetization
    Parameters
    ----------
    off
        offsets
    time
        time of the simulation
    B1
        B1 values to simulation
    use_pc
        to simulate phse cycling
    magn_init
        initial magnetization
    """
    if time is None:
        time = self.total_time
    # TODO code method/insert somwhere else
    # loop over phase cycling
    #   loop over offsets
    #       loop over B1
    #           loop over pulses
    #   phase cycling avg
    # return magn