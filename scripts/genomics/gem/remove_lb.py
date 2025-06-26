#!/usr/bin/env python

from pathlib import Path
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model, sbml
import cobra
import argparse

def clean_lb_reactions(input_path, output_path, remove=True):
    # Load the model
    model = cobra.io.read_sbml_model(input_path)

    # Process reactions
    for rxn in list(model.reactions):
        rxn_id_lower = rxn.id.lower()
        rxn_name_lower = (rxn.name or "").lower()
        if 'lb' in rxn_id_lower or 'lb' in rxn_name_lower:
            if remove:
                model.reactions.remove(rxn)
                print(f"Removed reaction: {rxn.id}")
            else:
                rxn.lower_bound = 0
                print(f"Blocked reaction: {rxn.id}")

    # Save the cleaned model, with error handling
    try:
        cobra.io.write_sbml_model(model, output_path)
        print(f"Saved cleaned model to {output_path}")
    except Exception as e:
        print(f"Error saving model: {e}")

def main():
    parser = argparse.ArgumentParser(description='Remove or block LB reactions in a SBML model.')
    parser.add_argument('input_file', help='Input SBML model file')
    parser.add_argument('output_file', help='Output SBML filename')
    parser.add_argument('--remove', action='store_true', help='Remove LB reactions entirely. Default: block them.')
    args = parser.parse_args()

    # Call the cleaning function
    clean_lb_reactions(args.input_file, args.output_file, args.remove)

if __name__ == '__main__':
    main()
