#in this one we try to implement bessel filter without scipy helpers

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

        

def set_poly_cutoff(poly, omega):
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
    print("")
    print("")
    n = 7
    omega1 = 30.0
    omega = np.logspace(0, 3, 100, base = 10)

    #generate a polynomial for bessel filter
    unscaled_poly_bessel = gen_poly(bessel_coef, n)
    print(f"n = {n}, unscaled poly = " + print_poly(unscaled_poly_bessel))
    
    #find its default 3dB frequency
    unscaled_omega_3db = get_poly_3db(unscaled_poly_bessel)
    print(f"n = {n}, default 3dB frequency = {unscaled_omega_3db:.8f} rad/s")
    
    #normalize to 3dB at 1
    poly_bessel = set_poly_cutoff(unscaled_poly_bessel, unscaled_omega_3db)
    print(f"n = {n}, scaled poly = " + print_poly(poly_bessel))
    
    #get a polynomial for butterworth filter
    poly_butter = gen_butter_poly(n)
    
    #move the cutoff of both to 10rads
    poly_bessel = set_poly_cutoff(poly_bessel, 1/omega1)
    poly_butter = set_poly_cutoff(poly_butter, 1/omega1)
    
    h0bes = get_poly_transfer(unscaled_poly_bessel, omega)
    h1bes = get_poly_transfer(poly_bessel, omega)
    h0but = get_poly_transfer(poly_butter, omega)
    
    gd0bes = get_gd_from_transfer(h0bes, omega)
    gd1bes = get_gd_from_transfer(h1bes, omega)
    gd0but = get_gd_from_transfer(h0but, omega)
    
    fig, (ax_mag, ax_phase, ax_gd) = plt.subplots(3, 1, figsize=(8, 6))
    plot_response(omega, h0bes, gd0bes, ax_mag, ax_phase, ax_gd, label="Unscaled Bessel", dotted=True)
    plot_response(omega, h1bes, gd1bes, ax_mag, ax_phase, ax_gd, label="Scaled Bessel")
    plot_response(omega, h0but, gd0but, ax_mag, ax_phase, ax_gd, label="Butterworth")
    ax_mag.axhline(-3, color='gray', linestyle='--', linewidth=0.7)

    plt.tight_layout()
    plt.legend()
    plt.show()
    
    
if __name__ == "__main__":
    main()