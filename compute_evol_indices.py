import os
import sys
import json
import argparse
import traceback
import pandas as pd
import torch

from EVE import VAE_model
from utils import data_utils

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Evol indices')
    parser.add_argument('--MSA_data_folder', type=str, help='Folder where MSAs are stored')
    parser.add_argument('--MSA_list', type=str, help='List of proteins and corresponding MSA file name')
    parser.add_argument('--protein_index', type=int, help='Row index of protein in input mapping file')
    parser.add_argument('--MSA_weights_location', type=str, help='Location where weights for each sequence in the MSA will be stored')
    parser.add_argument('--theta_reweighting', type=float, help='Parameters for MSA sequence re-weighting')
    parser.add_argument('--force_load_weights', action='store_true',
        help="Force loading of weights from MSA_weights_location (useful if you want to make sure you're using precalculated weights). Will fail if weight file doesn't exist.")
    parser.add_argument('--no_filter_columns', action='store_true', help='Do not filter columns by gap fraction')
    parser.add_argument('--VAE_checkpoint_location', type=str, help='Location where VAE model checkpoints will be stored')
    parser.add_argument('--model_name_suffix', default='Jan1', type=str, help='model checkpoint name is the protein name followed by this suffix')
    parser.add_argument('--model_parameters_location', type=str, help='Location of VAE model parameters')
    parser.add_argument('--computation_mode', type=str, help='Computes evol indices for all single AA mutations or for a passed in list of mutations (singles or multiples) [all_singles,input_mutations_list]')
    parser.add_argument('--all_singles_mutations_folder', type=str, help='Location for the list of generated single AA mutations')
    parser.add_argument('--mutations_location', type=str, help='Location of all mutations to compute the evol indices for')
    parser.add_argument('--output_evol_indices_location', type=str, help='Output location of computed evol indices')
    parser.add_argument('--output_evol_indices_filename_suffix', default='', type=str, help='(Optional) Suffix to be added to output filename')
    parser.add_argument('--num_samples_compute_evol_indices', type=int, help='Num of samples to approximate delta elbo when computing evol indices')
    parser.add_argument('--batch_size', default=256, type=int, help='Batch size when computing evol indices')
    parser.add_argument('--save_inprogress', action='store_true', help='Save and restore from in-progress batch checkpoints')
    args = parser.parse_args()

    mapping_file = pd.read_csv(args.MSA_list)
    protein_name = mapping_file['protein_name'][args.protein_index]
    msa_location = args.MSA_data_folder + os.sep + mapping_file['msa_location'][args.protein_index]
    print("Protein name: "+str(protein_name))
    print("MSA file: "+str(msa_location))

    if args.theta_reweighting is not None:
        theta = args.theta_reweighting
    else:
        try:
            theta = float(mapping_file['theta'][args.protein_index])
        except:
            theta = 0.2
    print("Theta MSA re-weighting: "+str(theta))

    if "weights_location" in mapping_file.columns:
            weights_location = args.MSA_weights_location + os.sep + mapping_file["weights_location"][args.protein_index]
            print("Using weights filename from mapping file")
    else:
        weights_location = args.MSA_weights_location + os.sep + protein_name + '_theta_' + str(theta) + '.npy'
    print("Weights file: "+str(weights_location))

    if args.force_load_weights:
        print("Flag force_load_weights enabled - Forcing that we use weights from file:", weights_location)
        if not os.path.isfile(weights_location):
            raise FileNotFoundError(f"Weights file {weights_location} doesn't exist."
                                    f"To recompute weights, remove the flag --force_load_weights.")

    data = data_utils.MSA_processing(
            MSA_location=msa_location,
            theta=theta,
            use_weights=True,
            weights_location=weights_location,
            threshold_focus_cols_frac_gaps=1.0 if args.no_filter_columns else 0.3,
    )
    
    if args.computation_mode=="all_singles":
        data.save_all_singles(output_filename=args.all_singles_mutations_folder + os.sep + protein_name + "_all_singles.csv")
        args.mutations_location = args.all_singles_mutations_folder + os.sep + protein_name + "_all_singles.csv"
    else:
        args.mutations_location = args.mutations_location + os.sep + protein_name + ".csv"
        
    model_name = protein_name + "_" + args.model_name_suffix
    print("Model name: "+str(model_name))

    model_params = json.load(open(args.model_parameters_location))

    model = VAE_model.VAE_model(
                    model_name=model_name,
                    data=data,
                    encoder_parameters=model_params["encoder_parameters"],
                    decoder_parameters=model_params["decoder_parameters"],
                    random_seed=42
    )
    model = model.to(model.device)

    checkpoint_name = str(args.VAE_checkpoint_location) + os.sep + model_name + "_final"
    try:
        checkpoint = torch.load(checkpoint_name)
        model.load_state_dict(checkpoint['model_state_dict'])
        print("Initialized VAE with checkpoint '{}' ".format(checkpoint_name))
    except Exception as e:
        print("Unable to restore from checkpoint '{}'".format(checkpoint_name))
        traceback.print_exc()
        sys.exit(1)

    output_filename = protein_name+'_'+str(args.num_samples_compute_evol_indices)+'_samples'+args.output_evol_indices_filename_suffix
    evol_indices_output_filename = args.output_evol_indices_location + os.sep + output_filename + ".csv"
    save_inprogress_prefix = None
    if args.save_inprogress:
        save_inprogress_location = "_inprogress"
        os.makedirs(save_inprogress_location, exist_ok=True)
        save_inprogress_prefix = save_inprogress_location + os.sep + output_filename

    list_valid_mutations, evol_indices, _, _ = model.compute_evol_indices(
        msa_data=data,
        list_mutations_location=args.mutations_location,
        num_samples=args.num_samples_compute_evol_indices,
        batch_size=args.batch_size,
        save_inprogress=args.save_inprogress,
        save_inprogress_prefix=save_inprogress_prefix,
    )

    df = {}
    df['protein_name'] = protein_name
    df['mutations'] = list_valid_mutations
    df['evol_indices'] = evol_indices
    df = pd.DataFrame(df)

    try:
        keep_header = os.stat(evol_indices_output_filename).st_size == 0
    except:
        keep_header=True 
    df.to_csv(path_or_buf=evol_indices_output_filename, index=False, mode='a', header=keep_header)
