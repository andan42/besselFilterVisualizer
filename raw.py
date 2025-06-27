#in this one we try to implement bessel filter without scipy helpers

from time import sleep

import matplotlib.pyplot as plt
import numpy as np


def factorial(n):
    if n == 0 or n == 1:
        return 1
    else:
        return n * factorial(n - 1)

def bessel_coef(k, n):
#sum of Ak*s^k
#where Ak =     fact(2n-k) / 
#               2^(n-k) * fact(k) * fact(n-k)
    if k > n:
        return 0
    else:
        return factorial(2 * n - k) / (2 ** (n - k) * factorial(k) * factorial(n - k))

def gen_butter_poly(n): #made by chatgpt replace this later with manual shit
    # Get complex conjugate pole pairs on left half-plane unit circle
    poles = [np.exp((-1j * np.pi * ((k + 1/2) / n + 1/2))) for k in range(n)]
    poly = np.poly(poles)  # converts roots → poly coefficients
    return poly.real.tolist()  # Should be real due to symmetry

def gen_poly(poly, n):
    return [poly(k, n) for k in range(n, -1, -1)]

def print_poly(poly):
    str_out = ""
    n = len(poly) - 1
    for i, coef in enumerate(poly[:-1]):
        if coef == 0:
            continue
        str_out += f"{coef:.2f} * x**{n - i} + "
    str_out += f"{poly[-1]:.2f}"
    return str_out

        

def scale_s_axis(poly, omega):
    n = len(poly) - 1
    return [coef * (omega ** (n - i)) for i, coef in enumerate(poly)]

def get_poly_3db(den, recursion_depth=1000, omega_start=0, omega_end=10000):
    omega = np.linspace(omega_start, omega_end, 1000)
    #avoid negative frequencies by replacing with zero
    omega = np.maximum(omega, 0)
    
    #find the index where the magnitude is closest to -3 dB
    h = den[-1] / np.polyval(den, 1j * omega)
    absh = np.abs(h)
    decibels = 20 * np.log10(absh)
    idx = np.argmin(np.abs(decibels + 3))
    #idx = np.argmin(np.abs(20 * np.log10(np.abs(h)) + 3))
    #if recursion is zero
    if recursion_depth <= 0 or omega[0] == omega[-1]:
        closest_omega = omega[idx]
        return omega[idx]
    elif idx == 0:
        next_omega_start = omega_start - 1
        next_omega_end = omega_end - 1
        return get_poly_3db(den, recursion_depth - 1, next_omega_start - 1, next_omega_end - 1)
    elif idx == len(omega) - 1:
        return get_poly_3db(den, recursion_depth - 1, omega_start + 1, omega_end + 1)
    else:
        closest_omega = omega[idx]
        current_delta_omega = omega_end - omega_start
        next_omega_start = closest_omega - current_delta_omega / 20
        next_omega_end = closest_omega + current_delta_omega / 20
        return get_poly_3db(den, recursion_depth - 10, next_omega_start, next_omega_end)

def get_poly_roots(poly):
    roots = np.roots(poly)
    return roots
    
def is_conjugate_pair(pole1, pole2):
    return np.iscomplex(pole1) and np.iscomplex(pole2) and np.isclose(np.conj(pole1), pole2)

def make_ordered_conjugate_pair(pole1, pole2):
    return (pole1, pole2) if np.imag(pole1) >= 0 else (pole2, pole1)

def is_pair_same(pairA, pairB):
    return (np.isclose(pairA[0], pairB[0]) and np.isclose(pairA[1], pairB[1])) or \
           (np.isclose(pairA[0], pairB[1]) and np.isclose(pairA[1], pairB[0]))
    
def pair_poles(poles, tol=1e-6):
    poles = list(poles)  # Ensure list, not np.array
    used = [False] * len(poles)
    paired = []
    unpaired = []

    for i, p1 in enumerate(poles):
        if used[i]:
            continue
        if np.iscomplex(p1):
            found_pair = False
            for j in range(i + 1, len(poles)):
                p2 = poles[j]
                if not used[j] and np.iscomplex(p2) and np.abs(p1 - np.conj(p2)) < tol:
                    paired.append([p1, p2])
                    used[i] = used[j] = True
                    found_pair = True
                    break
            if not found_pair:
                unpaired.append(p1)
                used[i] = True
        else:
            unpaired.append(p1)
            used[i] = True

    return paired, unpaired

def turn_pole_pair_to_poly(pole_pair):
    pole1, pole2 = pole_pair
    a2 = 1
    a1 = -(pole1 + pole2)
    a0 = pole1 * pole2
    return [a2, a1, a0]
        
def get_poly_transfer(den, omega):
    return den[-1] / np.polyval(den, 1j*omega)

def get_gd_from_transfer(h, omega):
    return -np.gradient(np.unwrap(np.angle(h))) / np.gradient(omega)

def plot_response(omega, h, gd, ax_mag, ax_phase, ax_gd, label=None, dotted=False):
    linestyle = '--' if dotted else '-'
    ax_mag.semilogx(omega, 20 * np.log10(np.abs(h)), label=label, linestyle=linestyle)
    ax_mag.set_ylabel("Magnitude [dB]")
    ax_mag.set_xlabel("Angular Frequency [rad/s]")
    ax_mag.grid(True)

    ax_phase.semilogx(omega, np.unwrap(np.angle(h)), label=label, linestyle=linestyle)
    ax_phase.set_ylabel("Phase [rad]")
    ax_phase.set_xlabel("Angular Frequency [rad/s]")
    ax_phase.grid(True)
    
    ax_gd.semilogx(omega, get_gd_from_transfer(h, omega), label=label, linestyle=linestyle)
    ax_gd.set_ylabel("Group Delay [s]")
    ax_gd.set_xlabel("Angular Frequency [rad/s]")
    ax_gd.grid(True)
    
    #plt.tight_layout()
    #plt.show()
   
def main():

    #Params for ECG filter
    n = 5 #Order. 2 MFB/SKs and a RC stage give 5 poles
    cutoff_filter = 2 * np.pi * 150  # 10 Hz in rad/s. Set 3dB point here

    #Response space
    #omega = np.logspace(-1, 2, 100, base=10)

    poly_bess_unnorm = gen_poly(bessel_coef, n)
    cutoff_bess_unnorm_den = get_poly_3db(poly_bess_unnorm)
    poly_bess_norm_den = scale_s_axis(poly_bess_unnorm, cutoff_bess_unnorm_den)
    poly_bess_den = scale_s_axis(poly_bess_norm_den, 1 / cutoff_filter)
    poles_bess = get_poly_roots(poly_bess_den)
    pair_poles_bess, unpair_pole_bess = pair_poles(poles_bess)
    mfb_poly_bess = [turn_pole_pair_to_poly(pair) for pair in pair_poles_bess]
  
    #organised print for stuff to turn into filters physically
    print("")
    print("")
    sleep(0.5)
 
    print("\nUnpaired poles:")
    for pole in unpair_pole_bess:
        print(f"Pole: {pole:.6f}")
    print("\nPaired poles:")
    for poly in mfb_poly_bess:
        print(f"Poly: {print_poly(poly)}")
    

    
    
    #plot the poles on a complex plane
    # plt.figure(figsize=(6, 6))
    # plt.scatter(np.real(poles_bess), np.imag(poles_bess), color='blue', label='Bessel Poles')
    # plt.axhline(0, color='black', linewidth=0.5, linestyle='--')
    # plt.axvline(0, color='black', linewidth=0.5, linestyle='--')
    
    # plt.grid(True)
    # plt.legend()
    # plt.show()
    
    

def test_transfers():
    print("")
    print("")
    n = 7
    omega1 = 7.0
    omega = np.logspace(-1, 2, 100, base = 10)

    #generate a polynomial for bessel filter
    unscaled_poly_bessel = gen_poly(bessel_coef, n)
    print(f"n = {n}, unscaled poly = " + print_poly(unscaled_poly_bessel))
    
    #find its default 3dB frequency
    unscaled_omega_3db = get_poly_3db(unscaled_poly_bessel)
    print(f"n = {n}, default 3dB frequency = {unscaled_omega_3db:.8f} rad/s")
    
    #normalize to 3dB at 1
    poly_bessel = scale_s_axis(unscaled_poly_bessel, unscaled_omega_3db)
    print(f"n = {n}, scaled poly = " + print_poly(poly_bessel))
    
    #get a polynomial for butterworth filter
    poly_butter = gen_butter_poly(n)
    
    #move the cutoff of both to 10rads
    poly_bessel = scale_s_axis(poly_bessel, 1/omega1)
    poly_butter = scale_s_axis(poly_butter, 1/omega1)
    
    #get resposes
    h0bes = get_poly_transfer(unscaled_poly_bessel, omega)
    h1bes = get_poly_transfer(poly_bessel, omega)
    h0but = get_poly_transfer(poly_butter, omega)
    
    #TODO but useless get their cutoffs.
    
    #get group delays
    gd0bes = get_gd_from_transfer(h0bes, omega)
    gd1bes = get_gd_from_transfer(h1bes, omega)
    gd0but = get_gd_from_transfer(h0but, omega)
    
    fig, (ax_mag, ax_phase, ax_gd) = plt.subplots(3, 1, figsize=(8, 6))
    plot_response(omega, h0bes, gd0bes, ax_mag, ax_phase, ax_gd, label="Unscaled Bessel", dotted=True)
    plot_response(omega, h1bes, gd1bes, ax_mag, ax_phase, ax_gd, label="Scaled Bessel")
    plot_response(omega, h0but, gd0but, ax_mag, ax_phase, ax_gd, label="Butterworth")
    ax_mag.axhline(-3, color='gray', linestyle='--', linewidth=0.7)
    ax_mag.axvline(omega1, color='red', linestyle='--', linewidth=0.7, label=f"Cutoff at {omega1:.2f} rad/s")
    ax_phase.axvline(omega1, color='red', linestyle='--', linewidth=0.7)
    ax_gd.axvline(omega1, color='red', linestyle='--', linewidth=0.7)

    plt.tight_layout()
    plt.legend()
    plt.show()
    
    
if __name__ == "__main__":
    main()