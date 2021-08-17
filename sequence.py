""" This module contains the class Sequence."""

import numpy as np
import pulse
from magnetization import simulate, magn_phase, plot_magn
import matplotlib.pyplot as plt


def pc_rec(phase_cycles, CTP):
    """
    Computes receiver phase with the phase cycling pc
    Parameters
    ----------
    phase_cyles: numpy array of floats
        array containing the phase cyle of each pulse on its line
    CTP: 1D-numpy array of floats
        vector containing the coherence transfer pathway of the pulses
    TODO code function
    """


class Sequence:

    """
    Class representing a pulse sequence
    """

    def __init__(self, pulses, total_time: float = None, pc=None):
        """
        Parameters
        ----------
        pulses: list of pulse objects
            the pulses of the sequence
        total_time: float
            Total duration of the sequence (s)
        pc: numpy array of floats
            phase cycling tables of the sequence
        """
        self.pulses = pulses

        if total_time is None:
            self.total_time = pulses[-1].end
        else:
            if total_time >= pulses[-1].end:
                self.total_time = total_time
            else:
                raise ValueError('total_time should be longer than the last'
                                 ' pulse end.')

        self.w1max = 0
        for p in pulses:
            if self.w1max < p.w1:
                self.w1max = p.w1

        if pc is not None:
            self.pc = pc

        # TODO: dicuss addition of list of pulse duration/delays
        # (can be deduced from start/stop in other functions)

    def __add__(self, sequence2add):
        """
        TODO code method
        """
        # check overlap?
        # pulses = [self.pulses, sequence2add.pulses]
        # sequence_sum = sequence(pulses, total_time=total_time)
        # return sequence_sum
        
    def __eq__(self, seq):
        """
        Parameters
        ----------
        seq: Sequence object
            sequence to compare
        Returns
        True/False: boolean            
        
        2 sequences are considered equal if they have the same total duration,
        the same pulses and all of their pulses are placed at the same 
        position. 
        """
                
        if self.total_time != seq.total_time:
            return False
       
        if self.pulses != seq.pulses:
            return False
        
        for p1, p2 in zip(self.pulses, seq.pulses):
            if p1.start != p2.start:
                return False
        
        return True

    def insert(self, elem, position):
        """
        Parameters
        ----------
        elem: pulse object or float
            element to insert in the sequence (pulse or delay)
        position: int
            where the pulse/delay should be inserted in the pulse list
        TODO code method
        """
        if isinstance(elem, pulse.Pulse):
            # insert pulse in pulse list
            # update total time and pulse positions
            pass
        elif isinstance(elem, float):  # delay case
            # update total time and pulse positions
            pass
        else:
            raise ValueError('A pulse (Pulse objects) or a delay (float) '
                             'should be inserted.')
            
    def append(self, elem):
        """
        Parameters
        ----------
        elem: pulse object or float
            element to insert in the sequence (pulse or delay)
        """
        self.insert(elem, len(self.pulses))

    def pulses_sum(self):
        """
        Returns the sum of the sequence pulses
        """
        return sum(self.pulses)

    def plot(self):
        """
        Plots the sequence pulses
        """
        p = sum(self.pulses)
        # TODO add potential final delay
        p.plot()

    def seq2TopSpin(self):
        """
        Export sequence pulse shapes to TopSpin waveforms
        TODO code method
        """

    def seq2Xepr(self):
        """
        Export sequence pulse shapes to Xepr waveforms
        TODO code method
        """
        # export each pulse with phase cycling options
        # plot each pulse and export its figure


class Exc_3fs(Sequence):

    """
    Class representing an excitation from 3 frequency-swept pulses (CHORUS,
                                                                   ABSTRUSE).
    # class or function?
    """

    def __init__(self, t90min, t180min, bw, tres, 
                 Q_exc=0.4412712003053, Q_ref=5, pulse_args={}, t_del=0,
                 polyfit=False, polyfit_args={}, 
                 plot=False):

        p1 = pulse.Parametrized(tp=t90min, bw=bw, Q=Q_exc, 
                                tres=tres, **pulse_args)
        p2 = pulse.Parametrized(tp=t180min+t90min/2, bw=bw, Q=Q_ref, 
                                tres=tres, start=t90min, **pulse_args)
        p3 = pulse.Parametrized(tp=t180min, bw=bw, Q=Q_ref, 
                                tres=tres, start=p2.end + t90min/2 + t_del,
                                **pulse_args)

        pulses = [p1, p2, p3]

        pc1 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        pc2 = np.array([0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3])
        pc3 = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3])

        #pc31 = pc_rec([pc1, pc2, pc3], [-1, +2, -2])

        pc31 = np.array([0, 2, 0, 2, 2, 0, 2, 0, 0, 2, 0, 2, 2, 0, 2, 0])
        pc = np.pi/2*np.stack((pc1, pc2, pc3, pc31))

        # simulation
        if polyfit:
            limit = 0.5*pulses[0].bw
            off = np.linspace(-limit, limit, 51)
            magn, off = simulate(pulses, offsets=off, pc=pc)

            if plot: plot_magn(magn, off)

            p1.add_ph_polyfit(magn_phase(magn), **polyfit_args, plot=plot)

        Sequence.__init__(self, pulses, p3.end+t_del, pc)
        # add other parameters?

        if plot:
            # simulation
            limit = 0.5*pulses[0].bw
            off = np.linspace(-limit, limit, 100)
            magn, off = simulate(self.pulses, offsets=off, pc=self.pc)
            
            plot_magn(magn, off)
            
            plt.figure()
            self.plot()
