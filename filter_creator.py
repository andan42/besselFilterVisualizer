from components_2 import (
    get_capacitor_values,
    get_my_capacitor_box_values,
    get_my_resistor_box_values,
    get_resistor_values,
    sort_and_filter_mfb_designs,
)
from raw import create_analog_filter, print_poly


def main():
    #use key value
    filter_params_list = {
        "Bessel 5th Order 150hz": 
        {
            "cutoff_hz": 150.0,
            "filter_order": 5,
            "filter_type": "bessel",
            "title": "Bessel 5th Order 150hz",
        },
        "Butterworth 5th Order 150hz": 
        {
            "cutoff_hz": 150.0,
            "filter_order": 5,
            "filter_type": "butterworth",
            "title": "Butterworth 5th Order 150hz",
        },
    }
    
    for filter_params in filter_params_list.values():
        do_it_all(
            cutoff_hz=filter_params["cutoff_hz"],
            filter_order=filter_params["filter_order"],
            filter_type=filter_params["filter_type"],
            title=filter_params["title"],
            gain_interval=(1.0, 2.0),  # Adjust as needed
        )
        
    return

def do_it_all(cutoff_hz, filter_order, filter_type, title, gain_interval=(None, None)):
    # cutoff_hz = 150.0 #in hz
    # filter_order = 5
    # filter_type = "bessel"
    #resistor_values = get_resistor_values([1e3, 1e4, 1e5])  # 1kΩ to 100kΩ
    #capacitor_values = get_capacitor_values([1e-9, 1e-8, 1e-7])  # 1nF to 100nF
    
    resistor_values = get_my_resistor_box_values()  # Use my actual resistor box values
    capacitor_values = get_my_capacitor_box_values()  # Use my actual capacitor box values
    
    print()
    print()
    print()
    print ("Filter Polynomials and Values for " + title)
    filter_poly_den_list = create_analog_filter(
        n_order = filter_order,
        cutoff_hz = cutoff_hz,
        filter_type = filter_type,
        plot_response_flag = False,
        plot_poles_flag = False,
        print_polys_flag = True,
    )

    for poly_den in filter_poly_den_list:
        if len(poly_den) != 3:
            continue
        a1_target = poly_den[1]
        a0_target = poly_den[2]
        df = sort_and_filter_mfb_designs(a1_target, 
                                    a0_target, 
                                    resistor_values=resistor_values, 
                                    capacitor_values=capacitor_values, 
                                    max_results=5, 
                                    gain_interval=gain_interval)
        
        # Print
        print("Component values for poly " + print_poly(poly_den) + " in filter: " + title)
        print()
        print(df)
        print()

    
    return

if __name__ == "__main__":
    main()