import itertools
from typing import List

import numpy as np

# Define E12 series standard values (scaled for base values)
E12_series = np.array([1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2])
decades = [1e3, 1e4, 1e5]  # 1kΩ to 100kΩ for resistors
capacitor_decades = [1e-9, 1e-8, 1e-7]  # 1nF to 100nF for capacitors

# E12_series = np.array([1.0, 2.2, 3.3, 4.7, 8.2])
# decades = [1e3, 1e4]  # 1kΩ to 100kΩ for resistors
# capacitor_decades = [1e-9, 1e-8]  # 1nF to 100nF for capacitors

# Build component value lists
resistor_values = np.sort(np.concatenate([E12_series * d for d in decades]))
capacitor_values = np.sort(np.concatenate([E12_series * d for d in capacitor_decades]))

def compute_a_terms(R1, R2, R3, C1, C2):
    a1 = (1 / C1) * (1/R1 + 1/R2 + 1/R3)
    a0 = 1 / (R1 * R2 * C1 * C2)
    return a1, a0

def compute_a_error(a1_target, a0_target, a1_eff, a0_eff):
    return np.sqrt(((a1_eff - a1_target) / a1_target)**2 + ((a0_eff - a0_target) / a0_target)**2)

#todo make compute gain function

def find_best_mfb_matches(a1_target, a0_target, max_results=5):
    results = []
    
    for R1, C1 in itertools.product(resistor_values, capacitor_values):
        for R2 in resistor_values:
            for C2 in capacitor_values:
                # Derive R3 from a1 equation (symbolically solved):
                a1_resid = a1_target - (1 / C1) * (1/R1 + 1/R2)
                if a1_resid <= 0:
                    continue
                R3 = 1 / (C1 * a1_resid)
                R3_rounded = min(resistor_values, key=lambda x: abs(x - R3))
                
                # Compute effective a1 and a0 with rounded R3
                a1_eff, a0_eff = compute_a_terms(R1, R2, R3_rounded, C1, C2)
                err = compute_a_error(a1_target, a0_target, a1_eff, a0_eff)
                
                results.append({
                    "R1": R1, "R2": R2, "R3": R3_rounded,
                    "C1": C1, "C2": C2,
                    "a1_eff": a1_eff, "a0_eff": a0_eff,
                    "gain (R3/R1)": R3_rounded / R1,
                    "error": err
                })

    # Sort by error and return top N
    sorted_results = sorted(results, key=lambda x: x["error"])[:max_results]
    return sorted_results

# Example placeholder target values (replace with user input in real use)
def main():
    example_a1 = 70000  # example units: 1/s
    example_a0 = 2.5e10  # example units: 1/s^2

    best_matches_df = find_best_mfb_matches(example_a1, example_a0)
    
    print("Best MFB Matches:")
    for match in best_matches_df:
        print(f"R1: {match['R1']:.2f} Ω, R2: {match['R2']:.2f} Ω, R3: {match['R3']:.2f} Ω, "
              f"C1: {match['C1']:.2e} F, C2: {match['C2']:.2e} F, "
              f"a1_eff: {match['a1_eff']:.2f}, a0_eff: {match['a0_eff']:.2e}, "
              f"Gain (R3/R1): {match['gain (R3/R1)']:.2f}, Error: {match['error']:.4f}")


if __name__ == "__main__":
    main()