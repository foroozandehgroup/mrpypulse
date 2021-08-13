
"""

This module contains functions for simulation purposes.

simulate uses rx0, ry0, rz0 and rtot0 imported from old JB project
simualte2 uses the functions developped by Ali (currently not working)

"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('seaborn-deep')


def simulate(pulses, offsets=None, tend=None, pc=None):
    """
    Calculation of the magnetisations Ix, Iy, and Iz along with the phase across a given 
    spectral width.
    
    Parameters
    ----------
    pulses: list of pulses
        pulses to be applied to magnetization
    offsets: numpy array of floats
        list of offsets at which to perform the simulation
    pc: numpy array of floats
        phase cycling to be used on pulses
        
    TODO: additional parameters
    nspins:
        if offsets=None
    B1: numpy array of floats
        B1 values to simulation
    magn_init: numpy array of floats
        initial magnetization
    TODO: single pulse case?
    
    Returns
    -------
    magn: ndarray
        ndarray containing the magnetisation across the spectral width        
    """
    if offsets is None:
        limit = 0.75*pulses[0].bw
        offsets = np.linspace(-limit, limit, 100)
        
    if tend is None:
        tend = pulses[-1].end
        
    if pc is None:
        pc = np.zeros((len(pulses)+1, 1))
    
    npc = pc.shape[1]
    noff = len(offsets)
    
    magn = np.zeros((3, len(offsets), npc))

    for phase in range(npc):
        for o in range(noff):
            
            # TODO: add option to start with different magentization vector
            
            # default: magnetization on z
            M = np.array([0, 0, 1])
            
            # potential delay before first pulse
            if pulses[0].start > 0:
                M = np.dot(rz0(2 * np.pi * offsets[o] * pulses[0].start), M)

            for i, p in enumerate(pulses):
                # if offsets[o]==0:
                #     print(pc[i,phase])
                
                # vectorized rotational matrix (contains all the pulse points)
                p_rtot = rtot0(p.r, offsets[o], pc[i,phase] + p.ph, p.tres)
                
                # apply rotational matrixes to magnetization
                for j in range(p.ns):                   
                    M = np.dot(p_rtot[:,:,j],M)
                    
                # potential delay between pulses
                if i < len(pulses) - 1:
                    
                    if p.end < pulses[i+1].start: # TODO np.isclose
                        M = np.dot(rz0(
                            2 * np.pi * offsets[o] * (pulses[i+1].start-p.end))
                            , M)

            # potential delay after last pulse
            if pulses[-1].end < tend:
                M = np.dot(rz0(
                    2 * np.pi * offsets[o] * (tend - pulses[-1].end)), M)

            magn[:, o, phase] = np.dot(rz0(-pc[-1,phase]), M)
            
    magn = np.sum(magn, axis=2) / npc # phase cycling collapse

    return magn, offsets


def simulate2(pulses, offsets=None, tend=None, pc=None):
    """
    Calculation of the magnetisations Ix, Iy, and Iz along with the phase across a given 
    spectral width.
    
    Parameters
    ----------
    pulses: list of pulses
        pulses to be applied to magnetization
    offsets: numpy array of floats
        list of offsets at which to perform the simulation
        
    # TODO: additional parameters
    nspins:
        if offsets=None
    B1: numpy array of floats
        B1 values to simulation
    pc: numpy array of floats
        phase cycling to be used on pulses
    magn_init: numpy array of floats
        initial magnetization
    
    Returns
    -------
    magn: ndarray
        ndarray containing the magnetisation across the spectral width  
        
    WARNING: currently not working
    """
    # Phase cycling
    # TODO: should be input by the user if he wants to use phase cycling


    if offsets is None:
        limit = 0.75*pulses[0].bw
        offsets = np.linspace(-limit, limit, 100)
        
    if tend is None:
        tend = pulses[-1].end
        
    if pc is None:
        pc = np.zeros((len(pulses)+1, 1))
    
    npc = pc.shape[1]
    noff = len(offsets)
    
    magn = np.zeros((3, len(offsets), npc))

    for phase in range(npc):
        for o in range(noff):
            
            # TODO: add option to start with different magentization vector
            M = np.array([0, 0, 1]) # default: magnetization on z
            
            # potential delay before first pulse
            if pulses[0].start > 0:
                M = np.dot(Rz(2 * np.pi * offsets[o] * pulses[0].start), M)

            for i, p in enumerate(pulses):
                
                for j in range(p.ns):
                    R_phi = pc[i,phase] + p.ph[j]
                    
                    # to use with 3*3 matrix
                    M = np.dot(Rtot(p.r[j], offsets[o], R_phi, p.tres), M)
                    
                    # to use with vector
                    # M = mvdot(Rtot(p.r[j], offsets[o], R_phi, p.tres), M)
                    
                # potential delay between pulses
                if i < len(pulses) - 1:
                    
                    if p.end < pulses[i+1].start: # TODO np.isclose
                        M = np.dot(Rz(2 * np.pi * offsets[o] * (pulses[i+1].start-p.end)), M)

            # potential delay after last pulse
            if pulses[-1].end < tend:
                M = np.dot(Rz(2 * np.pi * offsets[o] * (tend - pulses[-1].end)), M)

            magn[:, o, phase] = np.dot(Rz(-pc[-1,phase]), M)
            
    magn = np.sum(magn, axis=2) / npc # phase cycling collapse

    return magn, offsets


def magn_phase(magn):
    """
    Compute the phase of the magnetisation
    
    Parameters
    ----------
    magn: numpy array of floats
    
    Returns
    -------
    phase: 1D numpy array of floats
        phase of the magnetisation
    
    # TODO non 1D case
    """
    return np.angle(magn[1,:] + 1j * magn[0,:])


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


def plot_magn(magn, offsets):
    """
    Plot the magnetization.
    
    To prevent confusion of inputs and best results this should be called straight after the optimised pulse has been created,
    and from within the same script.
    TODO: change name? (some ppl use inversion/refocusing profiles...), we 
    could just rename to simulate()
    
    Parameters
    ----------
    magn: numpy array of floats
        magnetization to plot
    offsets: numpy array of floats    
        
    Plot over the desired bandwith the Ix, Iy, Iz, Ixy.
    
    """
        
    # # set figure parameters
    # mpl.rcParams['axes.spines.right'] = False
    # mpl.rcParams['axes.spines.top'] = False
    # mpl.rcParams['xtick.top'] = False
    # mpl.rcParams['ytick.right'] = False    

    # plot magnetisations against offset
    plt.figure() # 'Magnetisations vs offset'

    Ix = magn[0,:]
    Iy = magn[1,:]
    Iz = magn[2,:]
    Ixy = np.sqrt(Ix**2 + Iy**2)
    phase = np.angle(Iy + 1j * Ix)

    display = {"Ix": Ix, "Iy": Iy, "Iz": Iz, "Ixy":Ixy, "Phase":phase}

    nb_subplot = 1

    c = ['k', 'm', 'c', 'y', 'r']

    plt.suptitle('Magnetisations and phase vs offset')

    for [key, value], color in zip(display.items(), c):
        plt.subplot(5, 1, nb_subplot)
        plt.plot(offsets, value, color=color)
        plt.ylim(-1, 1)
        plt.ylabel(key)
        nb_subplot += 1
        if key == "Phase":
            plt.ylim(min(value), max(value))

    plt.xlabel("Offset")
    plt.tight_layout()

    return None 


def B1_variability(pulse, timestep, delta_F, max_rf_amp, nspins, N):
    '''
    This function plots the variability of the phase gainst B1 of the pulse. 

    To prevent confusion of inputs and best results this should be called straight after the optimised pulse has been created,
    and from within the same script. 

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
    '''
    #Plot phase variation with B1
    plt.figure('B1 phase variation')
    z, phi = B1_phase_variance(pulse, timestep, delta_F, max_rf_amp, nspins, N)

    plt.plot(z, phi, 'k')
    plt.xlabel(r'$\frac{B_1}{B_1^0}$')
    plt.ylabel(r'$\frac{d\phi}{dB_1}$')
    plt.title('Variation of phase across B1 frequencies')
    
    return None


def mvdot(matrix, vector):
    '''
    Matrix-vector multiplication rewritten for a column-flattened matrix
    and a row vector.

    Parameters
    ----------
    matrix: ndarray
        Column-flattened 3x3 matrix
    vector: ndarray
        Row vector

    Return
    ------
    result: ndarray
        Product of the multiplication as a row vector.
    '''
    
    if matrix.ndim == 1:
        matrix = matrix.reshape((1, -1))

    if vector.size == 3 and matrix.size == 9:
        result = np.zeros(3)
        result[0] = ((matrix[:,0] * vector[0]) + 
                    (matrix[:,3] * vector[1]) + 
                    (matrix[:,6] * vector[2]))

        result[1] = ((matrix[:,1] * vector[0]) + 
                    (matrix[:,4] * vector[1]) + 
                    (matrix[:,7] * vector[2]))

        result[2] = ((matrix[:,2] * vector[0]) + 
                    (matrix[:,5] * vector[1]) + 
                    (matrix[:,8] * vector[2]))

    else:
        result = np.zeros((np.size(vector, axis=0)))
        result[:,0] = ((matrix[:,0] * vector[:,0]) + 
                    (matrix[:,3] * vector[:,1]) + 
                    (matrix[:,6] * vector[:,2]))

        result[:,1] = ((matrix[:,1] * vector[:,0]) + 
                    (matrix[:,4] * vector[:,1]) + 
                    (matrix[:,7] * vector[:,2]))

        result[:,2] = ((matrix[:,2] * vector[:,0]) + 
                    (matrix[:,5] * vector[:,1]) + 
                    (matrix[:,8] * vector[:,2]))

    return result


def Rz(theta):
    '''
    Function to calculate the rotation by an angle theta around the z axis.
    The rotation matrix is returned as a column flattened vector.

    Parameters
    ----------
    theta: float or ndarray
        The angle(s) of rotation around the z axis

    Return
    ------
    Rz: ndarray
        numpy array of the rotation matrix

    '''
    cos = np.cos(theta)
    sin = np.sin(theta)

    rz = [cos, -sin, 0, sin, cos, 0, 0, 0, 1]
    
    return np.reshape(np.array([rz]), (3,3))


def Rtot(omega, offsets, phi, delta_time):
        '''
        Function to calculate the total rotation using the Rodrigues rotation formula. 
        Rotation around the effective magnetic field by the flip angle beta.
        Each rotation matrix is returned as a column flattened vector.

        Parameters
        ----------
        omega: float 
            RF amplitude
        offsets: float
            Offset frequency of the spin
        phi: float or ndarray
            Phase of pulse point
        delta_time: float
            Timestep

        Return
        ------
        rtot: ndarray
            numpy array of the rotation matrix
        '''

        if type(phi) == np.ndarray:
            omega_eff = np.zeros((len(omega), 2))
            omega_eff[:,0] = omega.reshape((-1,))
            omega_eff[:,1] = offsets
            omega_eff = np.linalg.norm(omega_eff * 2 * np.pi, axis=1) # effective field Beff

            theta = np.arctan2(omega, offsets) #Angle between B1 and Beff
            beta = delta_time * omega_eff #Flip angle
            beta = beta.reshape((-1, 1))

            phi_col = phi.reshape((-1,))
            N = len(phi)
            Beff = np.zeros((N,3))
            K = np.zeros((N, 9)) 
            K2 = np.zeros((N, 9))
            identity = np.eye(3).reshape((1, -1), order='F')
            Beff[:,0] = np.sin(theta) * np.cos(phi_col) 
            Beff[:,1] = np.sin(theta) * np.sin(phi_col)
            Beff[:,2] = np.cos(theta)

            Beff = Beff / np.linalg.norm(Beff, axis=1).reshape((-1,1))

            Beff2 = Beff ** 2 
            Beffxy = Beff[:,0] * Beff[:,1]
            Beffyz = Beff[:,1] * Beff[:,2]
            Beffzx = Beff[:,2] * Beff[:,0]

            K[:,0] = 0
            K[:,1] = Beff[:,2]
            K[:,2] = -Beff[:,1]
            K[:,3] = -Beff[:,2]
            K[:,4] = 0
            K[:,5] = Beff[:,0]
            K[:,6] = Beff[:,1]
            K[:,7] = -Beff[:,0]
            K[:,8] = 0
            K2[:,0] = -Beff2[:,2] - Beff2[:,1]
            K2[:,1] = Beffxy
            K2[:,2] = Beffzx
            K2[:,3] = Beffxy
            K2[:,4] = -Beff2[:,2] - Beff2[:,0]
            K2[:,5] = Beffyz
            K2[:,6] = Beffzx
            K2[:,7] = Beffyz
            K2[:,8] = -Beff2[:,1] - Beff2[:,0]

            rtot = identity + (np.sin(beta) * K) + ((1 - np.cos(beta)) * K2)
        
        else:

            omega_eff = np.sqrt((2 * np.pi * omega)**2       # Angular frequency of the
                        + (2 * np.pi * offsets)**2)  # effective field Beff
            theta = np.arctan2(omega, offsets) #Angle between B1 and Beff
            beta = delta_time * omega_eff #Flip angle
            
            Beff = np.zeros(3)
            K = np.zeros(9)
            K2 = np.zeros(9)
            identity = np.eye(3).reshape((1, -1), order='F')
            Beff[0] = np.sin(theta) * np.cos(phi) 
            Beff[1] = np.sin(theta) * np.sin(phi)
            Beff[2] = np.cos(theta)

            Beff = Beff / np.linalg.norm(Beff)

            Beff2 = Beff ** 2 
            Beffxy = Beff[0] * Beff[1]
            Beffyz = Beff[1] * Beff[2]
            Beffzx = Beff[2] * Beff[0]

            K[0] = 0
            K[1] = Beff[2]
            K[2] = -Beff[1]
            K[3] = -Beff[2]
            K[4] = 0
            K[5] = Beff[0]
            K[6] = Beff[1]
            K[7] = -Beff[0]
            K[8] = 0
            K2[0] = -Beff2[2] - Beff2[1]
            K2[1] = Beffxy
            K2[2] = Beffzx
            K2[3] = Beffxy
            K2[4] = -Beff2[2] - Beff2[0]
            K2[5] = Beffyz
            K2[6] = Beffzx
            K2[7] = Beffyz
            K2[8] = -Beff2[1] - Beff2[0]

            rtot = identity + (np.sin(beta) * K) + ((1 - np.cos(beta)) * K2)

        return rtot # np.transpose(np.reshape(rtot, (3,3)))
    

def rx0(phi):
    """Returns the rotational matrix for an angle phi around the x-axis
    
    If phi is an array containing n angles, return an array of n
    rotational matrixes for these angles around the x-axis.
    """
    if type(phi)==np.ndarray:

        m11 = np.full(len(phi),1)
        m12 = np.full(len(phi),0)
        m22 = np.cos(phi)
        m23 = np.sin(phi)

        m1 = np.stack((m11, m12, m12), axis=0)
        m2 = np.stack((m12, m22, m23), axis=0)
        m3 = np.stack((m12, m23, m22), axis=0)

        x_rot_mat = np.stack((m1, m2, m3), axis=0)

    else:
        x_rot_mat = np.array(([1, 0, 0], \
                            [0, np.cos(phi), np.sin(phi)], \
                            [0, np.sin(phi), np.cos(phi)]))

    return x_rot_mat


def ry0(phi):
    """Returns the rotational matrix for an angle phi around the y-axis
    """
    if type(phi)==np.ndarray:
        m11 = np.cos(phi)
        m12 = np.full(len(phi),0)
        m13 = np.sin(phi)
        m22 = np.full(len(phi),1)

        m1 = np.stack((m11, m12, m13), axis=0)
        m2 = np.stack((m12, m22, m12), axis=0)
        m3 = np.stack((-m13, m12, m11), axis=0)

        y_rot_mat = np.stack((m1, m2, m3), axis=0)

    else:
        y_rot_mat = np.array(([np.cos(phi), 0, np.sin(phi)], \
                            [0, 1, 0], \
                            [-np.sin(phi), 0, np.cos(phi)]))
    return y_rot_mat


def rz0(phi):
    """Returns the rotational matrix for an angle phi around the z-axis
    """
    if type(phi)==np.ndarray:
        m11 = np.cos(phi)
        m12 = np.sin(phi)
        m13 = np.full(len(phi),0)
        m33 = np.full(len(phi),1)

        m1 = np.stack((m11, -m12, m13), axis=0)
        m2 = np.stack((m12, m11, m13), axis=0)
        m3 = np.stack((m13, m13, m33), axis=0)

        z_rot_mat = np.stack((m1, m2, m3), axis=0)

    else:
        z_rot_mat = np.array(([np.cos(phi), -np.sin(phi), 0], \
                            [np.sin(phi), np.cos(phi), 0], \
                            [0, 0, 1]))
    
    return z_rot_mat


def rtot0(omega, offs, phi, time):
    """Returns the rotational matrix associated with a linear chirp
    inputs: pulse parameters
        - omega: point radiofrequency  or array containing points 
        associated with the B1 field
        - offs: point offset
        - phi: point phase
        - time
    output:
        - total_rot_mat: rotational matrix point associated with the 
        pulse point or array of rotational matrixes associated with the 
        pulse points
    Also computes omega_eff the angular frequency of the effective field
    Beff, theta the angle between Beff and B1 and the flip angle alpha 
    for the calculation of total_rot_mat
    """
    # Beff angular frequency
    omega_eff = np.sqrt((2 * np.pi * omega)**2 + (2 * np.pi * offs)**2)
    
    theta = np.arctan2(omega, offs) # angle between Beff and B1
    alpha = time * omega_eff # flip angle
    
    if type(phi)==np.ndarray:
        
        # array with each pulse point rotational matrix
        total_rot_mat = np.einsum('ijh,jkh,klh,lmh,mnh->inh',
            rz0(phi), ry0(theta), rz0(alpha), ry0(-theta),rz0(-phi))

    else:

        total_rot_mat = \
            rz0(phi).dot(ry0(theta).dot(rz0(alpha).dot(ry0(-theta).dot(rz0(-phi)))))
    
    return total_rot_mat
