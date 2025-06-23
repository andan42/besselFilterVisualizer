#%matplotlib inline  
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
from matplotlib.axes import Axes

mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams.update({'font.size': 12})


tau_accepted_dev = 0.05

def make_filter_bessel(
    order,
    f_cutoff,
    normalization = 'delay',
    n_points = 3000):

    w_cutoff = 2 * np.pi * f_cutoff


    [zeros, poles, gain] =  scipy.signal.bessel(order , w_cutoff, "lowpass", analog=True, output='zpk', norm=normalization)
    [num, den] =            scipy.signal.bessel(order , w_cutoff, "lowpass", analog=True, output='ba', norm=normalization)
    return num, den, zeros, poles, gain

def make_filter_butterworth(
    order,
    f_cutoff ,
    normalization = 'delay',
    n_points = 3000):
    
    w_cutoff = 2 * np.pi * f_cutoff
    [zeros, poles, gain] =  scipy.signal.butter(order , w_cutoff, "lowpass", analog=True, output='zpk')
    [num, den] =            scipy.signal.butter(order , w_cutoff, "lowpass", analog=True, output='ba')
    return num, den, zeros, poles, gain
    

def get_response(
    num,
    den,
    w_cutoff,
    n_points = 3000):
    #frequency response
    w = np.logspace(np.log10(w_cutoff / 1000), np.log10(w_cutoff * 10), 3000) 
    f = w / (2 * np.pi)

    w_resp, h = scipy.signal.freqs(num, den, worN = w)
    mag_db = 20 * np.log10(np.abs(h))
    phase = np.unwrap(np.angle(h))

    #group delay
    dphi = np.diff(phase)
    dw = np.diff(w_resp)
    tau = -dphi / dw
    f_tau = (f[:-1] + f[1:]) / 2
    
    num_d, den_d = scipy.signal.bilinear(num, den, fs=2000)
    
    _, tau_d = scipy.signal.group_delay(
        (num_d, den_d), 
        w=w_resp / (2 * np.pi), 
        fs=2000)
    tau_d = tau_d / 2000
    

    
    #todo do the stuff below
    # tau_dev = np.abs((tau - tau[0]) / tau[0])
    # idx = np.argmax(tau_dev > tau_accepted_dev)
    # f_flat_end = f_tau[idx] if idx else None
    return f, mag_db, phase, f_tau, tau, tau_d
def plot_poles_zeros(
    ax : Axes,
    zeros,
    poles,
    color
    ):
    
    ax.plot(np.real(zeros), np.imag(zeros), 'o', label='Zeros', color=color)
    ax.plot(np.real(poles), np.imag(poles), 'x', label='Poles', color=color)
    currentxlims = ax.get_xlim()
    currentylims = ax.get_ylim()
    leftmost = min(np.real(zeros).min() if len(zeros) !=0 else 0, 
                   np.real(poles).min() if len(poles) !=0 else 0)
    topmost = max(np.imag(zeros).max() if len(zeros) !=0 else 0,
                   np.imag(poles).max() if len(poles) !=0 else 0)
    bottommost = min(np.imag(zeros).min() if len(zeros) !=0 else 0,
                     np.imag(poles).min() if len(poles) !=0 else 0)
    if leftmost < currentxlims[0]:
        ax.set_xlim([leftmost, 0])
    if topmost > currentylims[1]:
        ax.set_ylim([bottommost, topmost])
    if bottommost < currentylims[0]:
        ax.set_ylim([bottommost, topmost])

def plot_response(
    ax,
    f,
    mag_db,
    phase,
    f_tau,
    tau,
    tau_d,
    color):
    
    ax[0].semilogx(f, mag_db, color=color) #, label=f'Order {len(tau)}'

    ax[1].semilogx(f, phase, color=color)

    ax[2].semilogx(f_tau, tau, label='τg', color=color, ls ='-', lw = 1)
    
    ax[2].semilogx(f, tau_d, color=color, ls='--', lw = 3)
    
    
    # ax[2].axhline(tau[0]*(1+tau_accepted_dev), color='grey', ls='--', lw=0.8)
    # ax[2].axhline(tau[0]*(1-tau_accepted_dev), color='grey', ls='--', lw=0.8)
    
    # if f_flat_end:
    #     ax[2].axvline(f_flat_end, color='red', ls='--', lw=0.8,
    #                   label=f'±{tau_accepted_dev*100:.0f}% τg limit ≈ {f_flat_end:.1f} Hz')

    # if f_flat_end: ax[2].legend()



def main():
    filters = [
        {'type': 'butterworth', 'order': 4, 'f_cutoff': 1/2*np.pi, 'normalization': 'mag', 'plot_color': 'orange'},
        {'type': 'bessel', 'order': 4, 'f_cutoff': 1/2*np.pi, 'normalization': 'mag', 'plot_color': 'blue'},
        
        #{'type': 'butterworth', 'order': 5, 'f_cutoff': 200, 'normalization': 'mag', 'plot_color': 'blue'},
        #{'type': 'butterworth', 'order': 5, 'f_cutoff': 300, 'normalization': 'mag', 'plot_color': 'red'},
        #{'type': 'butterworth', 'order': 5, 'f_cutoff': 400, 'normalization': 'mag', 'plot_color': 'green'},
    ]
         
    responses = []
    pole_zeros = []
    
    for m_filt in filters:
        if m_filt['type'] == 'bessel':
            num, den, zeros, poles, gain = make_filter_bessel(
                order=m_filt['order'], 
                f_cutoff=m_filt['f_cutoff'], 
                normalization=m_filt['normalization'])
        elif m_filt['type'] == 'butterworth':
            num, den, zeros, poles, gain = make_filter_butterworth(
                order=m_filt['order'], 
                f_cutoff=m_filt['f_cutoff'], 
                normalization=m_filt['normalization'])
        
        w_cutoff = 2 * np.pi * m_filt['f_cutoff']
        m_color = m_filt['plot_color']
        f, mag_db, phase, f_tau, tau, tau_d = get_response(num, den, w_cutoff)

        responses.append((f, mag_db, phase, f_tau, tau, tau_d, m_color))
        pole_zeros.append((zeros, poles, m_color))

    fig, ax = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
    
    fig2, ax2 = plt.subplots(1, 1)
    
    ax[0].set_ylabel('Mag (dB)')
    ax[0].grid(True, which='both')
    
    ax[1].set_ylabel('Phase (rad)')
    ax[1].grid(True, which='both')
    
    ax[2].set_ylabel('Group delay (s)')
    ax[2].grid(True, which='both')
    
    ax[2].set_xlabel('Frequency (Hz)')
    
    ax2.set_xlabel('Real')
    ax2.set_ylabel('Imaginary')
    ax2.set_title('Poles and Zeros')
    ax2.grid(False, which='both')
    

    #radius max is max of all poles and zeros magnitude divided by sqrt(2)
    max_complex_magnitude = 0
    for zeros, poles, _ in pole_zeros:
        if len(zeros) > 0:
            max_complex_magnitude = max(max_complex_magnitude, np.abs(zeros).max())
        if len(poles) > 0:
            max_complex_magnitude = max(max_complex_magnitude, np.abs(poles).max())
    radii = np.linspace(0, 2.0**(0.5) * max_complex_magnitude, 10)
    for r in radii:
        circ = mpatches.Circle((0, 0),    # centre at origin
                            radius=r,
                            fill=False,
                            lw=0.7, color='grey', alpha=0.6)
        ax2.add_patch(circ)
    angles = np.linspace(0, np.pi, 11)
    for angle in angles:
        ax2.plot([0, -np.sin(angle) * 2.0**(0.5) * max_complex_magnitude],
                 [0, np.cos(angle) * 2.0**(0.5) * max_complex_magnitude],
                 color='grey', lw=0.7, alpha=0.6)


    ax2.set_xlim([-1, 0])
    ax2.set_ylim([-1, 1])
    ax2.set_aspect(1, adjustable='box')
    
    for i, (f, mag_db, phase, f_tau, tau, tau_d, color) in enumerate(responses):
        plot_response(ax, f, mag_db, phase, f_tau, tau, tau_d, color)
        
    for i, (num, den, color) in enumerate(pole_zeros):
        plot_poles_zeros(ax2, num, den, color)
    
    #dont set any titles you god damn idiot
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.3) 
    plt.show()
    
if __name__ == "__main__":
    main()