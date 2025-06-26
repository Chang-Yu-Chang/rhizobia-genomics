#!/usr/bin/env python

from pathlib import Path
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model, sbml
import cobra

def set_tryptone_yeast_medium(model):
    # Close all exchange reactions initially
    # for reaction in model.exchanges:
    #     reaction.lower_bound = 0
    # Define key nutrients for Tryptone Yeast Medium (example metabolites)
    nutrients = {
        'EX_glc__D_e': -100,    # Glucose
        'EX_ala__L_e': -100,     # Alanine
        'EX_trp__L_e': -100,     # Tryptophan
        'EX_nad_e': -100,        # NAD+
        'EX_vitamin_b1_e': -100, # Vitamin B1
        # Add more metabolites as needed
    }
    for met_id, bound in nutrients.items():
        if met_id in model.reactions:
            model.reactions.get_by_id(met_id).lower_bound = bound
            #print(model.reactions.get_by_id(met_id).lower_bound)
        else:
            print(f"Reaction {met_id} not found in the model.")

for i in [2,3,4,5,15]:
    model = cobra.io.read_sbml_model('/Users/cychang2/Dropbox/lab/rhizobia-genomics/data/genomics/gem/xml/g'+str(i)+'.xml')  # update path
    # Output model information
    print("Model name: g" + str(i))
    print(f"Number of reactions: {len(model.reactions)}")
    print(f"Number of genes: {len(model.genes)}")
    # Set medium to Tryptone Yeast
    #set_tryptone_yeast_medium(model)
    # Run FBA
    solution = model.optimize()
    # Print the objective and fluxes
    print('Objective (biomass flux):', solution.objective_value)






import math

def thermal_performance_curve_C(T_celsius, T_opt_C=37, T_max_C=57, T_min_C=17, peak=1.0, drop_off=0.0):
    """
    Returns a performance factor (0 to 1) based on temperature in Celsius.
    - T_opt_C: Optimal temperature in °C
    - T_max_C: Max temperature in °C
    - T_min_C: Min temperature in °C
    """
    # Convert Celsius to Kelvin
    T = T_celsius + 273.15
    T_opt = T_opt_C + 273.15
    T_max = T_max_C + 273.15
    T_min = T_min_C + 273.15
    
    if T < T_min or T > T_max:
        return 0  # No activity outside bounds
    
    sigma = (T_max - T_min) / 6  # Controls the shape of the curve
    performance = peak * math.exp(-((T - T_opt) ** 2) / (2 * sigma ** 2))
    return max(performance, drop_off)

def apply_temp_sensitive_scaling_C(model, T_celsius):
    scale_factor = thermal_performance_curve_C(T_celsius)
    print(f"Applying scale factor {scale_factor:.2f} at {T_celsius}°C")
    for reaction in model.reactions:
        reaction.lower_bound *= scale_factor
        reaction.upper_bound *= scale_factor


for i in [2,3,4,5,15]:
    for t in [20, 30, 40]:
        model = cobra.io.read_sbml_model('/Users/cychang2/Dropbox/lab/rhizobia-genomics/data/genomics/gem/xml/g'+str(i)+'.xml')  # update path
        # Output model information
        print("\nModel name: g" + str(i) + " at ", str(t), "C")
        apply_temp_sensitive_scaling_C(model, t)
        # Run FBA
        solution = model.optimize()
        # Print the objective and fluxes
        print('Objective (biomass flux):', solution.objective_value)
