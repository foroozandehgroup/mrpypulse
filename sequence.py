""" This module contains the class Sequence."""

import numpy as np
import pulse
import simulation


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
            if total_time > pulses[-1].end:
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

    def insert(self, elem, position):
        """
        Parameters
        ----------
        elem: pulse object or float
            element to insert in the sequence (pulse or delay)
        position: int
            where the pulse/dealy should be inserted in the pulse list
        TODO code method
        """
        if isinstance(elem, pulse.Pulse):
            # insert pulse in pulse list
            # update total time and pulse positions
            pass
        if isinstance(elem, float):  # delay case
            # update total time and pulse positions
            pass

    def plot():
        """
        Plots the sequences pulses
        TODO code the method
        """
        # add all the pulse and plot the resulting total pulse?
        # limit: forced superposition but can be plotted by just ploting a
        # list of pulses on the same graph?

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

    def __init__(t90min, t180min, Q_exc=0.4412712003053, Q_ref=5, t_del=0,
                 plot=False, polyfit=False, pulse_param=(), polyfit_opt=(), **kwargs):

        p1 = pulse.Parametrized(t90min, Q_exc, pulse_param)
        p2 = pulse.Parametrized(1.25*t180min, Q_ref, pulse_param, start=t90min)
        p3 = pulse.Parametrized(t180min, Q_ref, pulse_param, start=p1.tp + p2.tp + t90min/2 + t_del)

        pulses = [p1, p2, p3]

        pc1 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        pc2 = np.array([0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3])
        pc3 = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3])

        pc31 = pc_rec([pc1, pc2, pc3], [-1, +2, -2])

        pc = [pc1, pc2, pc3, pc31]

        # simulation
        if polyfit:
            magn = simulation.magnetisations(pulses, p1.bw, 100)

            if plot:
                simulation.plot_magn(magn)

            p1.phase_polyfit(simulation.magn_phase(magn), polyfit_opt, plot)

        self = Sequence(pulses, p3.end+t_del, pc)
        # add other parameters?

        if plot:
            # simulation
            magn = simulation.magnetisations(pulses, p1.bw, 100)

            simulation.plot_magn(magn)

            self.plot()
