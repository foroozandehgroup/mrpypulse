""" This module contains the class Pulse and various inherited pulses."""

import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt

# TODO add Random()

class Pulse:

    """
    Class representing a pulse

    All attributes not input are generated.

    Required Parameters
    -------------------
        2 of tp, np and tres
        either x and y or r and ph

    Parameters
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
    ID: string
        Identification pulse
    TODO discuss addition of flip angle +others?
    TODO empty pulse case
    """

    def __init__(self, tp: float = None, ns: int = None, tres: float = None,
                 x=None, y=None,
                 r=None, ph=None,
                 phi0: float = 0, start: float = 0, ID: str = None):
        if ID is not None:
            self.ID = ID
        if x is None and y is None and r is None and ph is None:
            # no coordinates (empty pulse)
            pass
        elif x is None and y is None and r is not None and ph is not None:
            self.r = r
            self.ph = ph
            if len(r) != len(ph):
                raise ValueError('r and ph should have the same length')
        elif x is not None and y is not None and r is None and ph is None:
            self.x = x
            self.y = y
            if len(x) != len(y):
                raise ValueError('x and y should have the same length')
        else:
            raise TypeError('Coordinates should be input as'
                            ' either x and y or r and ph.')

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
            ns = tp/tres
            # self.ns = int(ns)
            if np.modf(ns)[0] > 0.99999:  # account for float type operations
                self.ns = int(np.ceil(ns))
            else:
                self.ns = int(np.floor(ns))
            self.tres = tp / self.ns  # tres adjusted for the rouding on ns
        else:
            raise TypeError('Exactly 2 of tp, np and tres should be used.')
        self.phi0 = phi0

        if self.ns > 10000:
            inputStr = input(
                'Number of pulse point > 10000! Input y to continue anyway.')
            if inputStr != "y":
                raise ValueError(
                    'High number of puls point, execution cancelled.')

        if hasattr(self, 'r'):

            if len(self.r) != self.ns:
                raise ValueError(
                    'Pulse coordinates length and ns do no match.')
            self.w1 = np.max(self.r)

        self.start = start
        self.end = self.start + self.tp

        # t defines the position of each segment (middle)
        self.t = np.linspace(
                 self.start+self.tres/2, self.end-self.tres/2, self.ns)

    def __setattr__(self, name, value):
        """
        Handles multiple attributes modification when one attribute is modified

        Set up the attribute identified by name with value. Attributes which
        causes other modifications:
        r or ph:
            x, y, w1 update
        x or y:
            r, ph, w1 update
        w1:
            scale coordinates x, y, and r
        phi0:
            phi0 also added to ph
        start, end :
            pulse position changed (start/end/t updated)
        t:
            needs to be compatible with pulse duration and number of points

        Particular conditions are set up to take into account the calls made in
        __init__
        """

        # already initialized special cases
        if name == 'w1' and hasattr(self, 'w1'):
            # scale coordinates (calls __setattr__ recursively)
            self.r = value * self.r/max(self.r)
        elif name == 'phi0' and hasattr(self, 'phi0'):
            # reset phi0 on ph
            object.__setattr__(self, 'ph', self.ph - self.phi0)

        elif name == 'start' and hasattr(self, 'start'):
            object.__setattr__(self, 'end', value + self.tp)
            object.__setattr__(
                self, 't',
                np.linspace(value+self.tres/2, self.end-self.tres/2, self.ns))

        elif name == 'end' and hasattr(self, 'end'):
            object.__setattr__(self, 'start', value - self.tp)
            object.__setattr__(
                self, 't',
                np.linspace(
                    self.start+self.tres/2, value-self.tres/2, self.ns))

        elif name == 't' and hasattr(self, 't'):
            object.__setattr__(self, 'start', value[0] - self.tres/2)
            object.__setattr__(self, 'end', value[-1] + self.tres/2)
            if np.isclose(value[-1]-value[0], self.tp, rtol=1e-6, atol=1e-15):
                raise ValueError('Vector t not compatible with tp.')
            if len(value) != self.ns:
                raise ValueError('t should has a different number of point.')

        # TODO tres/ns/tp modification and how they would affect each other
        # e.g.
        #   - tres -> modif np [requires interpolation]?
        #   - tp -> lower tres and no ns modif? +modify t/end
        #   - ns -> modif tres? [requires interpolation]?

        # set attribute value
        object.__setattr__(self, name, value)

        if name in ('ph', 'r') and hasattr(self, 'r') and hasattr(self, 'ph'):
            object.__setattr__(self, 'x', self.r * np.cos(self.ph))
            object.__setattr__(self, 'y', self.r * np.sin(self.ph))
            if hasattr(self, 'w1'):
                object.__setattr__(self, 'w1', np.max(self.r))

        elif name in ('x', 'y') and hasattr(self, 'x') and hasattr(self, 'y'):
            object.__setattr__(self, 'r', np.sqrt(self.x**2 + self.y**2))
            object.__setattr__(self, 'ph', np.arctan2(self.y, self.x))
            if hasattr(self, 'w1'):
                object.__setattr__(self, 'w1', np.max(self.r))

        elif name == 'phi0' and hasattr(self, 'phi0') and hasattr(self, 'ph'):
            self.ph = self.ph + self.phi0  # calls __setattr__ recursively

    def __add__(self, pulse2add):
        """
        Pulse addition operation defined as sum of Cartesian coordinates

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
        # TODO __rsub__
        
        if not np.isclose(self.tres, pulse2add.tres, rtol=1e-6, atol=1e-12):
            raise ValueError(
                'Pulses can only be added if their tres is the same.')

        # initialization
        tres = self.tres
        start = np.min([self.start, pulse2add.start])
        end = np.max([self.end, pulse2add.end])
        tp = end - start
        ns = tp/tres
        # TODO cleaner solution: only use int(), no rounding and catch errors if too far
        # not working tho
        if np.modf(ns)[0] > 0.99999:  # account for float type operations
            ns = int(np.ceil(ns))
        else:
            ns = int(np.floor(ns))

        ns = int(ns)
        tres = tp / ns  # tres adjusted for the rouding on n

        t = np.linspace(start+tres/2, end-tres/2, ns)

        x = np.empty(len(t))
        y = np.empty(len(t))
        j = 0  # first pulse index
        k = 0  # second pulse index

        for i in range(ns):

            x[i] = 0
            y[i] = 0

            if self.start < t[i] < self.end:
                x[i] += self.x[j]
                y[i] += self.y[j]
                j += 1

            if pulse2add.start < t[i] < pulse2add.end:
                x[i] += pulse2add.x[k]
                y[i] += pulse2add.y[k]
                k += 1

        pulse_sum = Pulse(
            ns=ns, tp=tp, x=np.array(x), y=np.array(y), start=start)

        return pulse_sum

    def __radd__(self, object2add):
        """
        Pulse add for non-pulses objects
        
        Parameters
        ----------
        object2add: object
            non-pulse object to add
        Returns
            self if the object is 0 (allows to use sum on an list of pulses)
        """
        if object2add == 0:
            return self
        else:
            raise ValueError('A pulse object should be added '
                             'to a pulse object.')

    def __eq__(self, p):
        """
        Parameters
        ----------
            - p, pulse to compare

        2 pulses are considered equal if they have the same coordinates,
        time resolution and number of points.
        They can take different positions and have different additional
        attributes.
        """
        eq = np.allclose(self.x, p.x, rtol=1e-6, atol=1e-15) and \
            np.allclose(self.y, p.y, rtol=1e-6, atol=1e-15) and \
            np.isclose(self.tres, p.tres, rtol=1e-6, atol=1e-15) and \
            self.ns == p.ns

        return eq

    def __ne__(self, p):
        """
        Negation of __eq__ (cf __eq__ for more information)
        """
        return not self.__eq__(p)

    def __str__(self):
        """
        Convert pulse object to string (typically used by print)
        """

        pulsestr = 'Pulse object with the following attributes\n'

        if hasattr(self, 'ID'):
            pulsestr += f'ID:    {self.ID}\n'

        pulsestr += (f'tp:      {self.tp}\n'
                     f'ns:      {self.ns}\n'
                     f'tres:    {self.tres}\n'
                     f'start:   {self.start}\n'
                     f'end:     {self.end}\n'
                     f'w1:      {self.w1}\n'
                     f'phi0:    {self.phi0}\n'
                     f'x:       [{self.x[0]} ... {self.x[-1]}]\n'
                     f'y:       [{self.y[0]} ... {self.y[-1]}]\n'
                     f'r:       [{self.r[0]} ... {self.r[-1]}]\n'
                     f'ph:      [{self.ph[0]} ... {self.ph[-1]}]\n')

        return pulsestr

    def add_ph_polyfit(self, ph, start=0, end=100, deg=5, plot=False):
        """
        Add the polynomial fitting of a phase vaector to the phase of the pulse

        Parameters
        ----------
        ph: numpy array of floats
            phase to add
        start: float
            start of the polynomial fit (%)
        end: float
            stop of the polynomial fit (%)
        deg: int
            degree of the polynomial fit
        plot: boolean
            allows to plot ph, its fit and the phase correction

        plt.show() might be needed calling the function to reveal the plots.

        """
        if start > end:
            raise ValueError('end should be superior to start')

        # select ph between start and end
        i_start = int(round(len(ph) * start / 100))
        i_end = int(round(len(ph) * end / 100))
        ph2fit = ph[i_start:i_end]

        # fit applied on selection
        x_ph2fit = np.arange(i_start, i_end)
        poly = np.polyfit(x_ph2fit, ph2fit, deg)
        fit = np.polyval(poly, x_ph2fit)

        # compute phase correction over whole pulse
        x_ph_corr = np.linspace(1, len(ph), self.ns)
        ph_corr = np.polyval(poly, x_ph_corr)     

        # apply phase correction
        self.ph += ph_corr

        if plot:
            plt.figure()
            plt.plot(x_ph2fit, ph2fit, x_ph2fit, fit, "r")

            plt.figure()
            plt.plot(x_ph_corr, ph_corr)

        return ph_corr

    def plot(self,
             form: str = "Cartesian", label: bool = True, title: str = None):
        """
        Plot the pulse shape in Cartesian coordinates

        Parameters
        ----------
        type: string
            type of plot
        title: string
            plot title

        Might require figure() before call and show() after.
        """
        if label:
            plt.xlabel('Time (s)')  # before potential call to twinx()

        if form == "Cartesian":
            plt.plot(self.t, self.x)
            plt.plot(self.t, self.y, 'r')

            plt.ylim(-self.w1, self.w1)
            if label:
                plt.ylabel('Cartesian coordinates (Hz)')

        elif form == "polar":
            plt.plot(self.t, self.r)
            plt.ylabel('Amplitude (Hz)', color='C0')
            plt.ylim(-self.w1, self.w1)

            ax = plt.gca()
            ax.twinx()
            plt.plot(self.t, self.ph, 'r')
            if label:
                plt.ylabel('Phase (rad)', color='r')
        else:
            raise ValueError('form should be one of the following: Cartesian, \
                             polar')

        plt.xlim(self.start, right=self.end)
        plt.title(title)

    def pulse2Xepr(self):
        """
        Export the pulse to Xepr format

        Returns
        -------
        x: numpy array of floats
            Pulse Cartesian coordinates x normalized from -1 to 1
        y: numpy array of floats
            Pulse Cartesian coordinates y normalized from -1 to 1
        """
        x_Xepr = self.x / self.w1
        y_Xepr = self.y / self.w1

        return x_Xepr, y_Xepr

    def pulse2TopSpin(self):
        """
        Export the pulse to TopSpin format

        Returns
        -------
        r_TopSpin: numpy array of floats
            Pulse amplitude normalized from 0 to 1
        ph_TopSpin: numpy array of floats
            Pulse phase normalized from 0° to 360°
        """
        r_TopSpin = self.r/self.w1
        ph_TopSpin = np.rad2deg(self.ph)

        return r_TopSpin, ph_TopSpin

    def pulse2Xepr_file(self):
        """
        Export the pulse to a shape file (.shp) for Xepr
        """

    def pulse2TopSpin_file(self):
        """
        Export the pulse to a shape file for TopSpin
        """

    def resonator_easyspin(self, eng, f, H_f, nu):
        """
        [EPR] Makes the pulse compensate for the resonator effects

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

        try:
            import matlab.engine
        except ImportError:
            print('matlab.engine could not be imported.')

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
        t2, y_t2 = eng.resonator(tmat, y_tmat, nu, fmat, H_fmat,
                                 'compensate', nargout=2)

        t2 = np.ravel(np.asarray(t2))
        y_t2 = np.ravel(np.asarray(y_t2))
        y_t2 = np.interp(t, t2, y_t2)

        y_t2 = np.real(y_t2) / max(np.real(y_t2)) + \
            1j * np.imag(y_t2) / max(np.imag(y_t2))

        self.r = 100 * abs(y_t2)
        self.ph = np.angle(y_t2, deg=True) % 360


class Hard(Pulse):

    """Class representing a hard pulse."""

    def __init__(self, tp, w1, **kwargs):
        """
        Parameters
        ----------
        tp : float
            cf. Pulse
        w1 : float
            cf. Pulse
        **kwargs
            other arguments to be transmitted to parent classes

        A hard pulse is defined as a 2 points pulse

        # TO DO: require testing
        """
        # TOD execution should raise an error
        Pulse.__init__(self, tp=tp, ns=2, tres=tp/2,
                       r=np.array([w1, w1]), ph=np.array([0, 0]),
                       **kwargs)


class Shape(Pulse):

    """
    Class representing a shaped pulse.

    Parameters
    ----------
    AM: string
        amplitude modulation type
    FM: string
        frequency modulation type
    bw: float
        bandwidth (Hz)
    **kwargs
        TODO document **kwargs throughout the code

    A shaped pulse is a pulse which can be amplitude-modulated (AM) and/or
    frequency-modulated
    """

    def __init__(self, AM: str = None, FM: str = None, bw: float = None,
                 **kwargs):

        Pulse.__init__(self, **kwargs)

        self.bw = bw
        
        if self.bw is not None:
            self.tbp = self.bw * self.tp
        else:
            self.tbp = None

        # test to distinguish no modulation from unknown modulation
        # (only if the pulse has coordinatess)
        if AM is None:           
            if hasattr(self, 'ph'):
                if np.all(self.r != self.r[0]):
                    self.AM = "unknown"
        if FM is None:
            if hasattr(self, 'ph'):
                if np.all(self.ph != self.ph[0]):
                    self.FM = "unknown"

    def __str__(self):
        """
        Convert shape object to string (typically used by print)
        """
        shape_str = super().__str__() +\
            (f'FM:      {self.FM}\n'
             f'AM:      {self.AM}\n'
             f'bw:      {self.bw}\n')

        return shape_str

    def reverse_sweep(self):
        """
        Reverse the sweep of a shape pulse

        Raises
        ------
        AttributeError
            if no FM is used, the pulse cannot be reversed
        """
        if self.FM is not None:
            self.y = -self.y
        else:
            raise AttributeError('No sweep to be reversed (FM=None).')


class Parametrized(Shape):

    """
    Class representing a parametrized AM/FM pulse

    Parameters
    ----------
    AM: string
        amplitude modulation, can take the following values: WURST,
        sinsmoothed (default), superGaussian, thanh
    FM: string
        frequency modulation, can take the following values: chirp (default),
        sech
    tp: float
        cf. Pulse
    bw: float
        cf. Shape
    w1: float
        cf. Pulse
    Q: float
        adiabaticity factor of the pulse
    delta_f: float
        frequency offset (by default 0, correspond to a centred FM)
    n: float
        smoothing index for WURST or superGaussian AM
    sm: float
        smoothing percentage for sinsmoothed AM
    B: float
        smoothing index for HS pulses
    **kwargs
        other argurments to transmit to parent classes

    Parametrized shaped pulses make use of analytical function for their
    waveforms.
    Exactly 3 of tp, bw, w1 and Q should be used to creat a parametrized shaped
    pulse which is frequency modulated.
    """

    def __init__(self, AM: str = "sinsmoothed", FM: str = "chirp",
                 tp: float = None, bw: float = None, w1: float = None,
                 Q: float = None, delta_f: float = 0,
                 n: int = None, sm: float = None, B: float = None, **kwargs):
        
        # required parameters
        if FM is not None:

            if w1 is None and bw is not None and \
                    tp is not None and Q is not None:
                w1 = np.sqrt(bw * Q / (2 * np.pi * tp))

            elif tp is None and bw is not None and \
                    Q is not None and w1 is not None:
                tp = bw * Q / (2 * np.pi * w1**2)

            elif bw is None and tp is not None and \
                    Q is not None and w1 is not None:
                bw = w1**2 * 2 * np.pi * tp / Q

            elif Q is None and bw is not None and \
                    tp is not None and w1 is not None:
                Q = w1**2 * 2 * np.pi * tp / bw

            else:
                raise TypeError('Exactly 3 of Q, w1, tp and bw should be used '
                                'as parameters for a frequency-modulated '
                                'pulse.')
            self.Q = Q
            self.w1 = w1
            
        elif AM is not None:
            if w1 is None:
                raise TypeError('w1 is needed for an amplitude-modulated '
                                'pulse.')
            else: self.w1=w1
        else:
            raise TypeError('No parametrized pulse with both AM and FM equal '
                            'to None')
            
        Shape.__init__(self, AM=AM, FM=FM, bw=bw, tp=tp, **kwargs)

        # frequency offset
        self.delta_f = delta_f

        # position delta_t
        self.delta_t = self.start + self.tp/2

        # amplitude modulation
        self.AM = AM
        if self.AM == "WURST":
            if n is None:
                n = 80  # default smoothing index value
            self.n = n

            self.r = self.w1 * (1 - np.abs(np.sin(
                     (np.pi * (self.t - self.delta_t)) / self.tp))**self.n)

        elif self.AM == "tanh":
            if B is None:
                B = 10.6/self.tp
            self.B = B

            # sech = 1/cosh
            self.r = self.w1 * 1/np.cosh(self.B*(self.t-self.delta_t))

        elif self.AM == "sinsmoothed":
            if sm is None:
                sm = 10  # default smoothing percentage value
            self.sm = sm

            # number of points smoothed
            n_sm = int(np.floor((self.ns * self.sm) / 100))

            # number of points unsmoothed
            n_unsm = int(self.ns - (2 * n_sm))

            unsmoothed_middle = self.w1 * np.ones(n_unsm)

            # amplitude apodized with a sine function taken from 0 to pi/2
            smoothed_side = self.w1 * (np.sin(np.linspace(0, np.pi/2, n_sm)))

            self.r = np.concatenate((smoothed_side,
                                     unsmoothed_middle,
                                     np.flip(smoothed_side)))

        elif self.AM == "superGaussian":
            if n is None:
                n = 40  # default smoothing factor (superGaussian index) n
            self.n = n
            self.r = self.w1 * \
                np.exp(
                    -2**(self.n + 2) *
                    ((self.t - self.delta_t) / self.tp)**self.n)

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
                      np.pi * self.bw * (self.t - self.delta_t)**2 / self.tp \
                      + 2 * np.pi * self.delta_f * (self.t - self.delta_t)

        elif self.FM == "sech":
            if B is None:
                B = 10.6/self.tp
            self.B = B

            # phase calculated from instantaneous frequency integral
            # instant_freq = 0.5 * bw * tanh(B*t);
            # instant_phase_integral = 0.5 * bw * (1/B) * log(cosh(B * (t)))
            self.ph = self.phi0 + \
                np.pi * self.bw * (1/self.B) * \
                np.log(np.cosh(self.B * (self.t - self.delta_t))) + \
                2 * np.pi * self.delta_f * (self.t - self.delta_t)

        elif self.FM is None:
            self.ph = self.phi0 * np.ones(self.ns)

    def __setattr__(self, name, value):
        """
        """
        # TODO delta_t modif
        # TODO other modifs with constructor?
        # TODO: tres/tp/ns -> modify them differntly
        # (no need to interpolate but would destroy pulse tweaks)

        # set attribute value
        Pulse.__setattr__(self, name, value)

    def __str__(self):
        """
        Convert parametrized object to string (typically used by print)
        """
        parametrized_str = super().__str__()

        if hasattr(self, 'FM'):
            parametrized_str += (f'delta_f: {self.delta_f}\n'
                                 f'Q:       {self.Q}\n')

        if hasattr(self, 'n'):
            parametrized_str += f'n:       {self.n}\n'
        if hasattr(self, 'sm'):
            parametrized_str += f'sm:      {self.sm}\n'
        if hasattr(self, 'B'):
            parametrized_str += f'B:       {self.B}\n'

        return parametrized_str

    def add_ph_polyfit(self, ph, start=None, end=None, deg=5, plot=False):
        """
        Parameters
        ----------
        """

        if start is None:
            if hasattr(self, 'sm'):
                start = self.sm
            else:
                start = 0

        if end is None:
            if hasattr(self, 'sm'):
                end = 100 - self.sm
            else:
                end = 100

        ph_corr = super().add_ph_polyfit(ph, 
                                      start=start, end=end, deg=deg, 
                                      plot=plot)

        return ph_corr


class OCT(Shape):

    """Class representing a Optimizal Control Theory (OCT) pulse - could be
    implemented in the future"""

    def __init__(self, algo=None, B1_range=None, **kwargs):
        Shape.__init__(self, **kwargs)


class Composite(Shape):

    """Class representing a composite pulse - could be implemented in the
    future"""

    pass
