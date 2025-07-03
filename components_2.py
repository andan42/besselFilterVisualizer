import bisect
import itertools

import numpy as np
import pandas as pd


def get_nearest_in_list(value, sorted_values):
    """
    Given a value and a sorted list, return a set of up to 2 nearest values:
    one below and one above. If value is out of bounds, return closest endpoint once.
    Output is a set, so duplicates are automatically removed.
    """
    index = bisect.bisect_left(sorted_values, value)

    if index == 0:
        # Below range
        return {sorted_values[0]}
    elif index >= len(sorted_values):
        # Above range
        return {sorted_values[-1]}
    else:
        # Between two values
        return {sorted_values[index - 1], sorted_values[index]}

def calculate_a_terms(R1, R2, R3, C1, C2):
    """Calculate a1, a0, and gain (G) from given component values."""
    a1 = (1 / C1) * (1 / R1 + 1 / R2 + 1 / R3)
    a0 = 1 / (R1 * R2 * C1 * C2)
    G = R3 / R1
    return a1, a0, G

def calculate_a_error(a1_target, a0_target, a1_eff, a0_eff):
    """Compute normalized error between target and effective a values."""
    rel_error = np.sqrt(((a1_eff - a1_target) / a1_target) ** 2 +
                        ((a0_eff - a0_target) / a0_target) ** 2)
    return rel_error

def solve_r1_r3(C1, C2, R2, a0_target, a1_target):
    """
    Given C1, C2, R2, a0 and a1, compute the ideal R1 and R3 values.
    Return None if the solution is invalid (e.g. R3 <= 0).
    """
    try:
        R1 = 1 / (a0_target * R2 * C1 * C2)
        inv_R3 = C1 * a1_target - (1 / R1 + 1 / R2)
        if inv_R3 <= 0:
            return None  # invalid R3
        R3 = 1 / inv_R3
        return R1, R3
    except ZeroDivisionError:
        return None

def generate_possible_mfb_designs(a1_target, a0_target, resistor_values, capacitor_values):
    """
here i will draw the MFB design in schematic form so you know which components are which. annoying cos theres 3 res and 2 caps so yeaaa

                                   R1
                    ---------------\/\/\-------------
                    |                               |
                    |                      C2       |
                    |               -------||-------|
                    |               |               |  
Vin ----\/\/\-------|------\/\/\-------POS          |
        R3          |      R2          OPAOPA       |
                    |                  OPAOPAOPA>---------------- Vout
                    |                  OPAOPA
                  ===== C1         ----NEG
                    |              |                              
                    |              |
                    =              =
                
                
idfk what else to do to make this clear whatever. 
input resistor is R3, input capacitor is C1
output feedback resistor is R1, output capacitor is C2
negative input reistor is R2
"""
    results = []

    # Iterate over C1, C2, and R2 combinations
    for C1, C2 in itertools.product(capacitor_values, capacitor_values):
        for R2 in resistor_values:
            # Solve ideal R1 and R3 for current C1, C2, R2
            solution = solve_r1_r3(C1, C2, R2, a0_target, a1_target)
            if solution is None:
                continue

            R1_ideal, R3_ideal = solution

            # Get sets of closest standard values
            R1_options = get_nearest_in_list(R1_ideal, resistor_values)
            R3_options = get_nearest_in_list(R3_ideal, resistor_values)

            # Try all valid combinations of rounded R1 and R3
            for R1_real in R1_options:
                for R3_real in R3_options:
                    #
                    if R1_real <= 0 or R3_real <= 0:
                        continue

                    # Compute effective a1, a0, and gain
                    a1_eff, a0_eff, G = calculate_a_terms(R1_real, R2, R3_real, C1, C2)
                    error = calculate_a_error(a1_target, a0_target, a1_eff, a0_eff)

                    results.append({
                        "C1": C1,
                        "C2": C2,
                        "R1": R1_real,
                        "R2": R2,
                        "R3": R3_real,
                        "Gain": G,
                        "a1_eff": a1_eff,
                        "a0_eff": a0_eff,
                        "Error": error
                    })

    # Sort results by error
    #results_sorted = sorted(results, key=lambda x: x["Error"])[:max_results]
    return results
    #return pd.DataFrame(results_sorted)

def sort_and_filter_mfb_designs(a1_target, a0_target, resistor_values, capacitor_values, max_results=10, gain_interval = (None, None)):
    """
    Generate and filter MFB designs based on target a1, a0, and gain.
    Returns a DataFrame of the best designs sorted by error.
    """
    results = generate_possible_mfb_designs(a1_target, a0_target, resistor_values, capacitor_values)

    # Convert to DataFrame for easier manipulation
    df = pd.DataFrame(results)

    # Filter by gain if specified
    if gain_interval[0] is not None:
        df = df[df['Gain'] >= gain_interval[0]]
    if gain_interval[1] is not None:
        df = df[df['Gain'] <= gain_interval[1]]

    # Sort by error and limit to max_results
    df = df.sort_values(by='Error').head(max_results)

    return df

def get_resistor_values(decades = [1e3, 1e4, 1e5]):
    """
    Generate a list of E12 series resistor values within a specified decade interval.
    Default is from 1kΩ to 100kΩ.
    """
    E12_series = np.array([1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2])
    return np.sort(np.concatenate([E12_series * d for d in decades]))
    
def get_my_resistor_box_values():
    """
    Hand typed box of actual resistors I have on hand.
    """
    return np.array([1.0e1, 2.2e1, 4.7e1,
                     1.0e2, 1.5e2, 2.0e2, 2.2e2, 2.7e2, 3.3e2, 4.7e2, 5.1e2, 6.8e2,
                     1.0e3, 2.0e3, 2.2e3, 3.3e3, 4.7e3, 5.1e3, 6.8e3,
                     1.0e4, 2.0e4, 4.7e4, 5.1e4, 6.8e4,
                     1.0e5, 2.2e5, 3.3e5, 4.7e5, 6.8e5,
                     1.0e6 ])
                     
def get_capacitor_values(decades = [1e-9, 1e-8, 1e-7]):
    """
    Generate a list of E12 series capacitor values within a specified decade interval.
    Default is from 1nF to 100nF.
    """
    E12_series = np.array([1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2])
    return np.sort(np.concatenate([E12_series * d for d in decades]))

def get_my_capacitor_box_values():
    """
    Hand typed box of actual capacitors I have on hand.
    """
    return np.array([1.0e-7, 1.5e-7, 2.2e-7, 3.3e-7, 4.7e-7, 6.8e-7,
                     1.0e-6, 2.2e-6, 4.7e-6, 
                     1.0e-5])
    
def all_pair_sums(available_values : np.ndarray):
    """
    Generate all possible sums that can be made by two added values from the available set.
    """
    value_pairs = itertools.combinations_with_replacement(available_values, 2)
    sum_values = [r1 + r2 for r1, r2 in value_pairs]
    return np.unique(sum_values)
    
def all_pair_inversesums(available_values : np.ndarray):
    """
    Generate all possible inverse sums that can be made by two added values from the available set.
    """
    value_pairs = itertools.combinations_with_replacement(available_values, 2)
    inversesum_values = [1 / (1/r1 + 1/r2) for r1, r2 in value_pairs]
    return np.unique(inversesum_values)

def all_pair_series_resistors(available_resistors : np.ndarray):
    """
    Generate all possible series combinations of resistors from the available set.
    """
    #use the helper
    return all_pair_sums(available_resistors)

def all_pair_parallel_resistors(available_resistors : np.ndarray):
    """
    Generate all possible parallel combinations of resistors from the available set.
    """
    #use the helper
    return all_pair_inversesums(available_resistors)

def all_pair_series_capacitors(available_capacitors : np.ndarray):
    """
    Generate all possible series combinations of capacitors from the available set.
    """
    #keep in mind correct series formula for capacitors is 1/C = 1/C1 + 1/C2
    return all_pair_inversesums(available_capacitors)

def all_pair_parallel_capacitors(available_capacitors : np.ndarray):
    """
    Generate all possible parallel combinations of capacitors from the available set.
    """
    #keep in mind correct parallel formula for capacitors is C = C1 + C2
    return all_pair_sums(available_capacitors)

def main():
    # Generate E12-based component value lists
    
    E12_series = np.array([1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2])

    # Define useful decades
    resistor_decades = [1e3, 1e4, 1e5]          # 1kΩ to 100kΩ
    capacitor_decades = [1e-9, 1e-8, 1e-7]      # 1nF to 100nF

    # Generate full value lists
    resistor_values = np.sort(np.concatenate([E12_series * d for d in resistor_decades]))
    capacitor_values = np.sort(np.concatenate([E12_series * d for d in capacitor_decades]))
    
    # Target values for a1 and a0
    
    a1_target = 2607.07 #1808.07  # example units: 1/s
    a0_target = 2158479.56 #2745846.14  # example units: 1/s
    
    # Search for valid MFB designs
    # results_df = generate_possible_mfb_designs(a1_target, a0_target, resistor_values, capacitor_values, max_results=10)
    results_df = sort_and_filter_mfb_designs(a1_target, a0_target, resistor_values, capacitor_values, max_results=10, gain_interval=(None, None))
    
    # Display results
    print("Found MFB designs:")
    print(results_df.to_string(index=False))
    
    # Compute A0 and A1 for the result table, as a sanity check
    print("\n")
    print("\n")
    print("\n")
    print("\nEffective a1, a0, and Gain for each design:")
    for index, row in results_df.iterrows():
        a1_eff, a0_eff, G = calculate_a_terms(row['R1'], row['R2'], row['R3'], row['C1'], row['C2'])
        print(f"Design {index + 1}: a1_eff = {a1_eff:.2f}, a0_eff = {a0_eff:.2f}, Gain = {G:.2f}, Error = {row['Error']:.6f}")
if __name__ == "__main__":
    main()