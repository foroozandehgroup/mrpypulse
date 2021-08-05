""" This module contains the class Sequence."""

class Sequence:
    
    "Class representing a pulse sequence"

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
                raise ValueError('total_time should be longer than the last pulse end.')

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
        """
        # TODO code method
        # check overlap?
        # pulses = [self.pulses, sequence2add.pulses]
        # sequence_sum = sequence(pulses, total_time=total_time)
        # return sequence_sum

    def insert_pulse(self, pulse, position):
        """
        Parameters
        ----------
        pulse: pulse object
            pulse to insert in the sequence
        position: int
            where the pulse should be inserted in the pulse list"""

        # TODO code method
        # insert pulse in pulse list
        # update total time and pulse positions

    def insert_delay(self, delay, position):
        """
        Parameters
        ----------
        delay: float
            delay to insert in the sequence
        position: int
            where the pulse should be inserted in the pulse list"""

        # TODO code method
        # update total time and pulse positions

    def pc_receiver(self):
        """Computes receiver phase with the phase cycling pc
        """
        # TODO code method/insert it somewhere else

    def seq2TopSpin(self):
        """Export sequence pulse shapes to TopSpin waveforms
        
        """
        # TODO code method

    def seq2Xepr(self):
        """Export sequence pulse shapes to Xepr waveforms"""
        # TODO code method
        # export each pulse with phase cycling options
        # plot each pulse and export its figure

