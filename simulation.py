
"This module contains functions for simulation purposes."

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('seaborn-deep')

# TODO this module

def simulate(pulses, off, time=None, B1=None, use_pc=False, magn_init=None):
    """
    Simulation of the magnetization
    
    Parameters
    ----------
    pulses: list of pulse objects
        pulses to be simulated
    off: numpy array of floats
        offsets
    time: float
        time of the simulation
    """
    if time is None:
        time = np.max(p.end for p in pulses)
    # TODO code method/insert somwhere else
    # loop over phase cycling
    #   loop over offsets
    #       loop over B1
    #           loop over pulses
    #   phase cycling avg
    # return magn


def magnetisations(pulses, delta_F, nspins):
    """
    Calculation of the magnetisations Ix, Iy, and Iz along with the phase across a given 
    spectral width.
    
    Parameters
    ----------
    pulses: list of pulses
        pulses to be applied to magnetization
    delta_F: float
        Spectral width in Hz
    nspins: int
        Number of spins
        
    # TODO: additional parameters
    B1: numpy array of floats
        B1 values to simulation
    pc: numpy array of floats
        phase cycling to be used on pulses
    magn_init: numpy array of floats
        initial magnetization
    
    Returns
    -------
    Ix: ndarray
        ndarray containing the Ix magnetisation across the spectral width
    Iy: ndarray
        ndarray containing the Iy magnetisation across the spectral width
    Iz: ndarray
        ndarray containing the Iz magnetisation across the spectral width
    phase: ndarray
        ndarray containing the phase across the spectral width
    offsets: ndarray
        ndarray containing offsets for which the values are calculated
        
    TODO use Rodrigues matrix for simulation
    """
    # Phase cycling
    # TODO: should be input by the user if he wants to use phase cyclign
    ph1 = [0]
    ph31 = [0]
    npc = len(ph1)

    # TODO: the user should input the offsets directly
    offsets = np.linspace(-0.75*delta_F, 0.75*delta_F, nspins)
    
    magnetisations = np.zeros((3, nspins, npc))

    for phase in range(npc):
        for o in range(len(offsets)):
            receiver = np.array([0, 0, 0])
            
            # TODO: add option to start differently
            M = np.array([0, 0, 1]) # Magnetisation starts as z-magnetisation
            
            # TODO loop over pulses and take delays into account
            for p in pulses:
                for i in range(p.ns):
                    #R_phi = np.pi * (ph1[phase] + p.ph[i]) / 180
                    R_phi = ph1[phase] + p.ph[i]
                    M = np.dot(Rtot(p.r[i], offsets[o], R_phi, p.tres), M)
            receiver = receiver + np.dot(Rz((-np.pi/180) * ph31[phase]), M)
            magnetisations[:, o, phase] = receiver

    Ix = magnetisations[0,:]
    Iy = magnetisations[1,:]
    Iz = magnetisations[2,:]
    Ixy = np.sqrt(Ix**2 + Iy**2)
    phase = np.angle(Iy + 1j * Ix) / np.pi

    # TODO: only return one magnetization array (other function for phase/Ixy)
    return Ix, Iy, Iz, Ixy, phase, offsets


def B1_phase_variance(pulse, timestep, delta_F, max_rf_amp, nspins, N):
    """
    Calculation of the magnetisations Ix, Iy, and Iz along with the phase across a given 
    spectral width.
    
    Parameters
    ----------
    Pulse: ndarray
        numpy array containing the control amplitudes (x_controls, y_controls). (1, 2N) array
    timestep: float
        Time of each pulse segment in seconds
    delta_F: float
        Spectral width in Hz
    max_rf_amp: float
        The power of the RF pulse in Hz
    nspins: int
        Number of spins
    N: int
        Number of pulse segments
    
    Returns
    -------
    z: ndarray
        ndarray containing the Iz magnetisation across the spectral width
    phi: ndarray
        ndarray containing the phase across the spectral width
        
    # TODO: adapt this function to pulse object
    """

    amp = np.sqrt((pulse[:N])**2 + (pulse[N:])**2)
    var = [0.995, 1.005]
    z = np.linspace(0.5, 1.5, nspins)

    #Phase cycling
    ph1 = [0]
    ph31 = [0]
    npc = len(ph1)

    rf_amps = max_rf_amp * amp

    phi = np.arctan2(pulse[N:], pulse[:N])
    phi = 180 + ((phi * 180) / np.pi)

    magnetisations = np.zeros((3, nspins, 2, npc))
    for spin in  range(nspins):
        for elem in range(2):
            for phase in range(npc):

                receiver = np.array([0, 0, 0])
                M = np.array([0,0, 1])
                for i in range(N):

                    R_phi = np.pi * (ph1[phase] + phi[i]) / 180
                    M = np.dot(Rtot(rf_amps[i]*var[elem]*z[spin], 0, R_phi, timestep), M)

                receiver = receiver + np.dot(Rz((-np.pi/180) * ph31[phase]), M)
                magnetisations[:, spin, elem, phase] = receiver

    Ix = magnetisations[0,:,:]
    Iy = magnetisations[1,:,:]
    
    phi = (180/np.pi) * (np.angle(Iy[:,1] + 1j * Ix[:,1]) 
                         - np.angle(Iy[:,0] + 1j * Ix[:,0]))

    for spin in range(nspins):
        if phi[spin] < -180:
            phi[spin] = phi[spin] + 360
        elif phi[spin] > 180:
            phi[spin] = phi[spin] - 360
        
    return z, phi


def excitation_profile(pulses, delta_F, nspins):
    """
    Calculates and plot the excitaton profile of the pulse.
    
    To prevent confusion of inputs and best results this should be called straight after the optimised pulse has been created,
    and from within the same script.
    TODO: change name? (some ppl use inversion/refocusing profiles...), we 
    could just rename to simulate()
    
    Parameters
    ----------
    Pulse: list of pulses
        pulses to simulate
    delta_F: float
        Spectral width in Hz
    nspins: int
        Number of spins
        
    Plot over the desired bandwith the Ix, Iy, Iz, Ixy.
    
    """
    
    #Set figure parameters
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False
    mpl.rcParams['xtick.top'] = False
    mpl.rcParams['ytick.right'] = False    

    #Plot magnetisations against offset
    plt.figure('Magnetisations vs offset')
    Ix, Iy, Iz, Ixy, phase, offsets = magnetisations(pulses, delta_F, nspins)

    display = {"Ix": Ix, "Iy": Iy, "Iz": Iz, "Ixy":Ixy, "Phase":phase}
    
    offsets = offsets * 1e-3
    nb_subplot = 1

    c = ['k', 'm', 'c', 'y', 'r']

    plt.suptitle('Magnetisations and phase vs offset')

    for [key, value], color in zip(display.items(), c):
        plt.subplot(5, 1, nb_subplot)
        plt.plot(offsets, value, color=color)
        plt.ylim(-1, 1)
        plt.ylabel(key)
        nb_subplot += 1

    plt.xlabel("Offset / kHz")
    plt.tight_layout()

    # TODO: move the follwing somewhere else (not an excitation profile)
    # Plot phase variation with B1
    # phase variation plotted vs offset
    # plt.figure('B1 phase variation')
    # z, phi = B1_phase_variance(pulse, timestep, delta_F, max_rf_amp, nspins, N)

    # plt.plot(z, phi, 'k')
    # plt.xlabel(r'$\frac{B_1}{B_1^0}$')
    # plt.ylabel(r'$\frac{d\phi}{dB_1}$')
    # plt.title('Variation of phase across B1 frequencies')

    # plt.show()
    return None # TODO: determine no return is better
    
    
def Rx(theta):
    """
    Function to calculate the rotation by an angle theta around the x axis
    
    Parameters
    ----------
    theta: float or ndarray
        The angle(s) of rotation around the x axis (rad)
    
    Returns
    ------
    Rx: ndarray
        numpy array of the rotation matrix
    """
    cos = np.cos(theta)
    sin = np.sin(theta)

    m1 = np.stack((1, 0, 0), axis=0)
    m2 = np.stack((0, cos, sin), axis=0)
    m3 = np.stack((0, sin, cos), axis=0)

    return np.stack((m1, m2, m3), axis=0)


def Ry(theta):
    """
    Function to calculate the rotation by an angle theta around the y axis
    Parameters
    ----------
    theta: float or ndarray
        angle(s) of rotation around the y axis (rad)
    
    Returns
    -------
    Ry: ndarray
        numpy array of the rotation matrix
    """
    cos = np.cos(theta)
    sin = np.sin(theta)

    m1 = np.stack((cos, 0, sin), axis=0)
    m2 = np.stack((0, 1, 0), axis=0)
    m3 = np.stack((-sin, 0, cos), axis=0)

    return np.stack((m1, m2, m3), axis=0)


def Rz(theta):
    """
    Function to calculate the rotation by an angle theta around the z axis
    
    Parameters
    ----------
    theta: float or ndarray
        angle(s) of rotation around the z axis (rad)
    
    Returns
    -------
    Rz: ndarray
        numpy array of the rotation matrix
    """
    cos = np.cos(theta)
    sin = np.sin(theta)

    m1 = np.stack((cos, -sin, 0), axis=0)
    m2 = np.stack((sin, cos, 0), axis=0)
    m3 = np.stack((0, 0, 1), axis=0)
    return np.stack((m1, m2, m3), axis=0)


def Rtot(omega, offsets, phi, delta_time):
    """
    Function to calculate the total rotation.
    
    Parameters
    ----------
    omega: float 
        amplitude (Hz)
    offsets: float
        offset frequency of the spin (Hz)
    phi: float or ndarray
        phase of pulse point (rad)
    delta_time: float
        timestep (s)
    
    Returns
    -------
    rtot: ndarray
        numpy array of the rotation matrix
    """

    omega_eff = np.sqrt((2 * np.pi * omega)**2       # Angular frequency of the
                        + (2 * np.pi * offsets)**2)  # effective field Beff

    theta = np.arctan2(omega, offsets) # Angle between B1 and Beff
    beta = delta_time * omega_eff      # Flip angle

    if type(phi)==np.ndarray:
        rtot = np.einsum('ijh,jkh,klh,lmh,mnh->inh',
                          Rz(phi), Ry(theta), Rz(beta), Ry(-theta), Rz(-phi))
    else:
        rtot = Rz(phi).dot(Ry(theta).dot(Rz(beta).dot(Ry(-theta).dot(Rz(-phi)))))
    
    return rtot

