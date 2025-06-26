import cobra
import argparse

def run_static_fba(model_path):
    # Load the model
    model = cobra.io.read_sbml_model(model_path)

    # Set media: example for glucose, adjust as needed
    # Initialize media
    concentrations = {
        'glc__D': 20,  # mmol/L
    }
    # Max uptake rates in mmol/gDW/hr
    max_uptake = {
        'glc__D': 10,
    }

    
    # For example, if you want to allow glucose uptake:
    rxn_id = 'EX_glc__D_e'  # change if your model has a different ID
    if rxn_id in [r.id for r in model.reactions]:
        model.reactions.get_by_id(rxn_id).lower_bound = -10  # allow uptake
        print(f"Set {rxn_id} lower_bound to -10")
    else:
        print(f"Reaction {rxn_id} not found in the model.")

    # Run FBA
    solution = model.optimize()

    # Output results
    print('Objective (biomass flux):', solution.objective_value)
    print('Flux distribution:')
    for r_id, flux in solution.fluxes.items():
        if abs(flux) > 1e-6:
            print(f"{r_id}: {flux:.4f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run static FBA on a model.')
    parser.add_argument('model_file', help='Path to your SBML model file.')
    args = parser.parse_args()

    run_static_fba(args.model_file)
