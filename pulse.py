""" This module contains the class Pulse and various inherited pulses."""

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
try:
    import matlab.engine
except:
    pass


class Pulse:

    """Class representing a pulse"""

    def __init__(self, tp=None, ns=None, tres=None, 
                 x=None, y=None, 
                 r=None, ph=None, 
                 phi0=0, start=0, ID=None):
        """
        All attributes not input are generated.

        Required Parameters
        -------------------
            2 of tp, np and tres
            either x and y or r and ph

        Attributes
        ----------
        tp: float
            duration (s)
        ns: int
            number of points/segments
        tres: float
            time resolution/sampling interval/segment duration
            (s/sample)
        x: numpy array of float
            first cartesian coordinates(Hz)
        y: numpy array of float
            second cartesian coordinates (Hz)
        r: numpy array of float
            amplitude polar coordinate (Hz)
        ph: numpy array of float
            phase polar coordinate (rad)
        phi0: float
            phase offset (rad)
        w1: float
            pulse maximum amplitude (Hz)
        t: numpy array of float
            an array containing the times of the pulse point (in s)
        start: float
            start of the pulse (s)
        end: float
            end of the pulse
        TODO discuss addition of flip angle and ID (name) +others?
        """
        if x is None and y is None and r is None and ph is None:
            # no coordinates (empty pulse)
            pass
        elif x is None and y is None:
            self.r = r
            self.ph = ph
            if len(r) != len(ph):
                raise ValueError('r and ph should have the same length')
        elif r is None and ph is None:
            self.x = x
            self.y = y
            if len(x) != len(y):
                raise ValueError('x and y should have the same length')
        else:
            raise ValueError('x should be input with y or r with ph.')

        if tres is None and tp is not None and ns is not None:
            self.tp = tp
            self.ns = ns
            self.tres = tp / ns
        elif tp is None and ns is not None and tres is not None:
            self.tp = ns * tres
            self.ns = ns
            self.tres = tres
        elif ns is None and tp is not None and tres is not None:
            self.tp = tp
            self.ns = int(np.floor(tp / tres))
            self.tres = tp / self.ns # tres adjusted for the rouding on n
        else:
            raise ValueError('Exactly 2 of tp, np and tres should be input.')

        self.phi0 = phi0

        if hasattr(self, 'r'):
            if len(self.r) != ns:
                raise ValueError('Pulse coordinates length and ns do no match.')
            self.w1 = np.max(self.r)

        self.start = start
        self.end = self.start + self.tp

        # t defines the position of each segment (defined in the middle)
        self.t = np.linspace(self.start+self.tres/2, self.end-self.tres/2, self.ns)

    def __setattr__(self, name, value):
        """
        Handles multiple attributes modification when one attribute is modified
        In particular, called in __init__
        """        
        if name == 'w1' and hasattr(self, 'w1'):
            # w1 already initialized: scale coordinates (calls __setattr__ recursively)
            self.r = value * self.r/max(self.r)     
        elif name == 'phi0' and hasattr(self, 'phi0'):
            # phi0 already initialized: reset phi0 on ph
            object.__setattr__(self, 'ph', self.ph - self.phi0)

        # TODO tres/ns/tp modification and how they would affect each other
        # e.g.
        #   - tres -> modif np [requires interpolation]?
        #   - tp -> lower tres and no ns modif?
        #   - ns -> modif tres?

        # set attribute value
        object.__setattr__(self, name, value)
        
        if (name == 'ph' or name == 'r') and hasattr(self, 'r') and hasattr(self, 'ph'):
            # r or ph: x, y, w1 update
            object.__setattr__(self, 'x', self.r * np.cos(self.ph))
            object.__setattr__(self, 'y', self.r * np.sin(self.ph))
            if hasattr(self, 'w1'):
                object.__setattr__(self, 'w1', np.max(self.r))
        elif (name == 'x' and name == 'y') and hasattr(self, 'x') and hasattr(self, 'y'):
            # x or y: r, ph, w1 update
            object.__setattr__(self, 'r', np.sqrt(self.x**2 + self.y**2))
            object.__setattr__(self, 'ph', np.arctan2(self.y, self.x))
            if hasattr(self, 'w1'):
                object.__setattr__(self, 'w1', np.max(self.r))
        elif name == 'phi0' and hasattr(self, 'phi0') and hasattr(self, 'ph'):
            # phi0: add phase offset to phase
            self.ph = self.ph + self.phi0 # calls __setattr__ recursively

    def __radd__(self, pulse2add):
        """Addition operation defined as sum of Cartesian coordinates
        
        Parameters
        ----------
        pulse2sum: pulse object
            pulse to sum

        Returns
        -------
        pulse_sum: pulse object
            sum of the pulses at their position

        If the pulses do not overlap, the delay between them is encoded as 
        coordinates of value 0.
        Attribute information from original pulses (such as phase offset)
        is lost.
        """
        # TODO testing
        if self.tres != pulse2add.tres:
            raise ValueError('Pulses can only be added if their tres is the same.')

        # initialization
        tres = self.tres
        start = np.min(self.start, pulse2add.start)
        end = np.max(self.end, pulse2add.end)
        tp = end - start
        t = np.arange(start+self.tres/2, end-self.tres/2, tres)
        x = np.empty(len(t))
        y = np.empty(len(t))
        j = 0 # first pulse index
        k = 0 # second pulse index

        for i, val in t:
            
            x[i]=0
            y[i]=0

            if self.start < val < self.stop:
                x[i] += self.x[j]
                j += 1

            if self.start < val < self.stop:
                x[i] += self.x[k]
                k += 1

        pulse_sum = Pulse(tres=tres, tp=tp, x=np.array(x), y=np.array(y), start=start)

        return pulse_sum

    def plot(self, name=""):
        """Plot the shape of the pulse in cartesian coordinates

        The figures needs to be revealed with plt.show() after the 
        function execution.

        """
        plt.figure()
        plt.plot(self.t, self.x)
        plt.plot(self.t, self.y, 'r')
        plt.ylim(-self.w1, self.w1)
        plt.xlim(self.start, right=self.end)
        plt.title(name)

    def pulse2Xepr(self):
        """Export the pulse to Xepr format"""

    def pulse2TopSpin(self):
        """Export the pulse to TopSpin format"""

    def pulse2Xepr_file(self):
        """Export the pulse to a shape file (.shp) for Xepr"""

    def pulse2TopSpin_file(self):
        """Export the pulse to a shape file for TopSpin"""
    
    def resonator_easyspin(self, eng, f, H_f, nu):
        """ [EPR] Makes the pulse compensate for the resonator effects
        
        It calls the resonator function from Easyspin with in Matlab
        Matlab, Easyspin and the Matlab engine are required:
        http://mathworks.com/products/matlab.html
        http://mathworks.com/help/matlab/matlab_external/matlab-engine-for-python.html
        http://www.easyspin.org
        For Matlab engine execution:
        eng = matlab.engine.start_matlab()
        pulse.resonator_easyspin(self, eng, f, H_f, nu)
        eng.quit()

        Parameters
        ----------
        eng: MatLabEngine object
            the MATLAB engine attached to Matlab process
        f: numpyt array of floats
            frequency vector of the transfer function H_f
        H_f: numpyt array of floats
            transfer function of the resonator
        nu: float
            the frequency on which to center the resonator compensation, in GHz
        
        The transfer function of the resonator H_f is used to change the
        pulse shape to account for the resonator.
        The pulse amplitude amp and phase ph are modified.
        """
        # TODO test this method
        t = self.t * 1e6
        y_t = self.x + self.y * 1j
        nres = f.size*8
        
        # Easyspin resonator function needs unique values
        H_f, index = np.unique(H_f, return_index=True) 
        f_interp = scipy.interpolate.interp1d(f[index], H_f, kind='cubic')

        # number of points of resonator profile extension
        f = np.linspace(min(f[index]), max(f[index]), nres)
        H_f = f_interp(f)

        # converting to Matlab variables
        tmat = matlab.double(t.tolist())
        y_tmat = matlab.double(y_t.tolist(), is_complex=True)
        fmat = matlab.double(f.tolist())
        H_fmat = matlab.double(H_f.tolist())

        # Easyspin resonator function call
        t2, y_t2 = eng.resonator(tmat, y_tmat, nu, fmat, H_fmat, \
                                 'compensate', nargout=2)

        t2 = np.ravel(np.asarray(t2))
        y_t2 = np.ravel(np.asarray(y_t2))
        y_t2 = np.interp(t, t2, y_t2)

        y_t2 = np.real(y_t2) / max(np.real(y_t2)) + 1j * np.imag(y_t2) / max(np.imag(y_t2))

        self.r = 100 * abs(y_t2)
        self.ph = np.angle(y_t2,deg=True) % 360


class Hard(Pulse):

    "Class representing a hard pulse."

    def __init__(self, tp, w1, **kwargs):
        Pulse.__init__(self, tp=tp, ns=2, tres=tp/2, 
                       r=np.array([w1, w1]), ph=np.array([0, 0]), 
                       phi0=0, start=0)


class Shape(Pulse):

    "Class representing a shaped pulse."

    def __init__(self, AM=None, FM=None, bw=None, sm=None, n=None, **kwargs):
       
        Pulse.__init__(self, **kwargs)

        self.bw = bw
        
        # test to distinguish no modulation from unknown modulation
        if AM is None:
            # test on self.r
            if np.all(self.r != self.r[0]):
                self.AM = "unknown"
        if FM is None:
            if np.all(self.ph != self.ph[0]):
                self.FM = "unknown"

    def reverse_sweep(self):
        if self.FM is not None:
            self.y = -self.y


class Parametrized(Shape):

    "Class representing a parametrized AM/FM pulse"

    def __init__(self, AM="sinsmoothed", FM="chirp", 
                 tp=None, bw=None, w1=None, Q=None, delta_f=0,
                 n=None, sm=None, B=None, **kwargs):
        """
        """

        if FM is not None:
            if w1 is None and bw is not None and tp is not None and Q is not None:
                w1 = np.sqrt(bw * Q / (2 * np.pi * tp))
            elif tp is None and bw is not None and Q is not None and w1 is not None:
                tp = bw * Q / (2 * np.pi * w1**2)
            elif bw is None and tp is not None and Q is not None and w1 is not None:
                bw = w1**2 * 2 * np.pi * tp / Q
            elif Q is None  and bw is not None and tp is not None and w1 is not None:
                Q = w1**2 * 2 * np.pi * tp / bw
            else:
                raise ValueError('Exactly 3 of Q, w1, tp and bw should be input.')
            self.Q = Q
            self.w1 = w1
            

        Shape.__init__(self, AM=AM, FM=FM, bw=bw, tp=tp, **kwargs)

        # frequency offset
        self.delta_f = delta_f

        # position delta_t
        self.delta_t = self.start + self.tp/2
        
        # amplitude modulation
        self.AM = AM
        if self.AM == "WURST":
            if n is None:
                n = 80 # default smoothing index value
            self.n = n

            self.r = self.w1 * (1 - np.abs(np.sin((np.pi * (self.t - self.delta_t)) / self.tp))**self.n)

        elif self.AM == "tanh":
            if B is None:
                B = 10.6/self.tp
            self.B = B

            # sech = 1/cosh
            self.r = self.w1 * 1/np.cosh(self.B*(self.t-self.delta_t))


        elif self.AM == "sinsmoothed":
            if sm is None:
                sm = 10 # default smoothing percentage value
            self.sm = sm

            # number of points smoothed
            n_sm = np.floor((self.ns * self.sm) / 100)

            # number of points unsmoothed
            n_unsm = self.ns - (2 * n_sm)

            unsmoothed_middle = self.w1 * np.ones(int(n_unsm))

            # amplitude apodized with a sine function taken from 0 to pi/2
            smoothed_side = self.w1 * (np.sin(np.linspace(0, np.pi/2, n_sm)))

            self.r = np.concatenate((smoothed_side, 
                                     unsmoothed_middle, 
                                     np.flip(smoothed_side)))

        elif self.AM == "superGaussian":
            if n is None:
                n = 40 # default smoothing factor (superGaussian index) n
            self.n = n

            self.r = self.w1 * np.exp(-(2**(self.n + 2)) * ((self.t - self.delta_t) / self.tp)**self.n)

        elif self.AM is None:
            self.r = self.w1 * np.ones(self.ns)

        """
        TODO
        # estimation of smoothing percentage sm if possible
        if self.AM == "WURST" or self.AM == "superGaussian"
            try:
                i_sm = 1 # unsmoothed part beginning index

                while self.r(i_sm) < 0.99 * p.w1
                    i_sm += 1
                end
                self.sm = 100 * i_sm / self.ns;
            catch:
                print('Smoothing sm could not be computed')
        """       

        # phase/frequency modulation
        self.FM = FM
        if self.FM is None:
            self.ph = np.zeros(self.ns)

        elif self.FM == "chirp":
            # phase calculated from instantaneous frequency
            # d(phase)/dt = sweep_rate * t + f0 (= instant. phase.)
            # sweep_rate = bw / tp;     
            # instant_phase_integral = (sweep_rate * t**2) / 2 + f0 * t;

            self.ph = self.phi0 + \
                      np.pi * self.bw * (self.t - self.delta_t)**2 / self.tp + \
                      2 * np.pi * self.delta_f * (self.t - self.delta_t)
        
        elif self.FM == "sech":
            if B is None:
                B = 10.6/self.tp
            self.B = B

            # phase calculated from instantaneous frequency integral
            # instant_freq = 0.5 * bw * tanh(B*t);
            # instant_phase_integral = 0.5 * bw * (1/B) * log(cosh(B * (t)))
            self.ph = self.phi0 + \
                      np.pi * self.bw * (1/self.B) * np.log(np.cosh(self.B * (self.t - self.delta_t))) + \
                      2 * np.pi * self.delta_f * (self.t - self.delta_t)
        
        elif self.FM is None:
            self.ph = self.phi0 * np.ones(self.ns)
    
    def __setattr__(self, name, value):
        """
        """
        # TODO delta_t modif
        # TODO other modifs with constructor?

        # set attribute value
        Pulse.__setattr__(self, name, value)

        
class OCT(Shape):

    "Class representing a Optimizal Control Theory (OCT) pulse - could be implemented in the future"

    def __init__(self, algo=None, B1_range=None, **kwargs):
        Shape.__init__(self, **kwargs)


class Composite(Shape):

    "Class representing a composite pulse - could be implemented in the future"

    pass

