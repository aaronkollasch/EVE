import os, sys
import glob
import argparse
import pandas as pd
import json
import torch

from EVE import VAE_model
from utils import data_utils

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='VAE')
    parser.add_argument('--MSA_data_folder', type=str, help='Folder where MSAs are stored')
    parser.add_argument('--MSA_list', type=str, help='List of proteins and corresponding MSA file name')
    parser.add_argument('--protein_index', type=int, help='Row index of protein in input mapping file')
    parser.add_argument('--MSA_weights_location', type=str, help='Location where weights for each sequence in the MSA will be stored')
    parser.add_argument('--theta_reweighting', type=float, help='Parameters for MSA sequence re-weighting')
    parser.add_argument('--no_filter_columns', action='store_true', help='Do not filter columns by gap fraction')
    parser.add_argument('--VAE_checkpoint_location', type=str, help='Location where VAE model checkpoints will be stored')
    parser.add_argument('--model_name_suffix', default='Jan1', type=str, help='model checkpoint name will be the protein name followed by this suffix')
    parser.add_argument('--model_parameters_location', type=str, help='Location of VAE model parameters')
    parser.add_argument('--training_logs_location', type=str, help='Location of VAE model parameters')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    parser.add_argument('--save_inprogress', action='store_true', help='Save and restore from in-progress model checkpoints')
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
        except (ValueError, KeyError):
            theta = 0.2
    print("Theta MSA re-weighting: "+str(theta))

    data = data_utils.MSA_processing(
            MSA_location=msa_location,
            theta=theta,
            use_weights=True,
            weights_location=args.MSA_weights_location + os.sep + protein_name + '_theta_' + str(theta) + '.npy',
            threshold_focus_cols_frac_gaps=1.0 if args.no_filter_columns else 0.3,
    )

    model_name = protein_name + "_" + args.model_name_suffix
    print("Model name: "+str(model_name))

    model_params = json.load(open(args.model_parameters_location))

    if args.save_inprogress:
        os.makedirs("_inprogress", exist_ok=True)
        model_params["training_parameters"]["save_model_params_freq"] = model_params["training_parameters"]["num_training_steps"] // 100
        model_params["training_parameters"]["model_inprogress_checkpoint_location"] = "_inprogress"
        old_checkpoints = glob.glob(model_params["training_parameters"]['model_inprogress_checkpoint_location']+os.sep+model_name+"_step_*")
        if len(old_checkpoints) > 0:
            args.seed += int(old_checkpoints[0].split("_")[-1])
    else:
        old_checkpoints = []

    model = VAE_model.VAE_model(
                    model_name=model_name,
                    data=data,
                    encoder_parameters=model_params["encoder_parameters"],
                    decoder_parameters=model_params["decoder_parameters"],
                    random_seed=args.seed
    )
    model = model.to(model.device)

    if len(old_checkpoints) > 0:
        checkpoint_name = max(old_checkpoints, key=lambda x: int(x.split("_")[-1]))
        try:
            checkpoint = torch.load(checkpoint_name)
            model.load_state_dict(checkpoint['model_state_dict'])
            print("Initialized VAE with checkpoint '{}' ".format(checkpoint_name))
        except:
            print("Unable to locate VAE model checkpoint")
            sys.exit(0)

    model_params["training_parameters"]['training_logs_location'] = args.training_logs_location
    model_params["training_parameters"]['model_checkpoint_location'] = args.VAE_checkpoint_location

    print("Starting to train model: " + model_name)
    model.train_model(data=data, training_parameters=model_params["training_parameters"])

    print("Saving model: " + model_name)
    model.save(model_checkpoint=model_params["training_parameters"]['model_checkpoint_location']+os.sep+model_name+"_final", 
                encoder_parameters=model_params["encoder_parameters"], 
                decoder_parameters=model_params["decoder_parameters"], 
                training_parameters=model_params["training_parameters"]
    )
