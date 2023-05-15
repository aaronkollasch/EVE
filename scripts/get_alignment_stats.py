#!/usr/bin/env python
import os
import subprocess
import pandas as pd
import sys
from functools import partial
from contextlib import contextmanager,redirect_stderr,redirect_stdout
sys.path.append("/home/ubuntu/EVE")
from utils import data_utils
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map

@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(os.devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)

def s3_cp_file(from_path, to_path, silent=False):
    max_runs = 3
    run = 0
    error = Exception()
    while run < max_runs:
        try:
            return subprocess.run(
                ["aws", "s3", "cp", from_path, to_path],
                check=True,
                stdout=subprocess.DEVNULL if silent else None,
                stderr=subprocess.DEVNULL if silent else None,
            )
        except subprocess.CalledProcessError as e:
            error = e
            continue
        else:
            break
        finally:
            run += 1
    else:
        raise error

s3_cp_file("s3://markslab-us-east-2/colabfold/output/popeve/missing_genes_with_priority.tsv", ".")
df = pd.read_table("missing_genes_with_priority.tsv")
assert len(df["protein"].unique()) == len(df)
prio_df = df[df['high_priority']]

os.makedirs("/tmp/aln_stats", exist_ok=True)

def get_aln_stats(run_name, dest_name, protein):
    s3_cp_file(f"s3://markslab-private/eve/popeve/data/MSA/{run_name}/{protein}.a2m", f"/tmp/aln_stats/{protein}.a2m")
    theta = 0.2
    weights_fname = protein + '_theta_' + str(theta) + '.npy'
    weights_location = "/tmp/aln_stats" + os.sep + weights_fname
    max_runs = 3
    run = 0
    error = Exception()
    weights_exist = False
    while run < max_runs:
        try:
            result = subprocess.run(
                [
                    "aws", "s3", "cp",
                    f"s3://markslab-private/eve/popeve/data/weights/{dest_name}/{weights_fname}",
                    weights_location,
                ],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
            )
        except subprocess.CalledProcessError as e:
            if b'does not exist' in e.stderr:
                weights_exist = False
                break
            else:
                error = e
                continue
        else:
            weights_exist = True
            break
        finally:
            run += 1
    else:
        raise error
    if not weights_exist:
        print(f"weights file s3://markslab-private/eve/popeve/data/weights/{dest_name}/{weights_fname} does not exist and will be created")
    with suppress_stdout_stderr():
        msa = data_utils.MSA_processing(
            f"/tmp/aln_stats/{protein}.a2m",
            weights_location=weights_location,
            threshold_focus_cols_frac_gaps=1.0,
        )
    assert msa.focus_seq_name.lstrip(">")  == protein
    num_seqs = 0
    with open(f"/tmp/aln_stats/{protein}.a2m") as f:
        for line in f:
            if line.startswith(">"):
                num_seqs += 1
    if not weights_exist:
        s3_cp_file(weights_location, f"s3://markslab-private/eve/popeve/data/weights/{dest_name}/")
    os.remove(f"/tmp/aln_stats/{protein}.a2m")
    os.remove(weights_location)
    result = {
        "protein": protein,
        "num_seqs": msa.num_sequences,
        "num_seqs_unfiltered": num_seqs,
        "perc_filtered": msa.num_sequences / num_seqs,
        "n_eff": msa.Neff,
        "n_eff_l": msa.Neff / len(msa.focus_seq),
        "seq_len": len(msa.focus_seq),
        "num_cov": msa.one_hot_encoding.shape[1],
        "perc_cov": msa.one_hot_encoding.shape[1] / len(msa.focus_seq),
    }
    del msa
    print(protein)
    return result

for run_name in [
#     "popeve_priority_proteins_uniref30_2103",
#     "popeve_priority_proteins_uniref30_2202",
#     "popeve_priority_proteins_uniref30_2202_c50",
#     "popeve_priority_proteins_uniref30_2202_c25",
    "popeve_priority_proteins_uniref30_2202_c40",
]:
    dest_name = run_name + "_m0"
    result = process_map(
        partial(get_aln_stats, run_name, dest_name),
        (row.protein for row in prio_df.itertuples()),
        max_workers=len(os.sched_getaffinity(0)),
        chunksize=10,
        total=len(prio_df),
    )
    result_df = pd.DataFrame(result)
    # result_df.to_csv(f"{run_name}_aln_stats.csv", index=False)
    result_df.to_csv(f"{dest_name}_aln_stats.csv", index=False)
    s3_cp_file(f"./{dest_name}_aln_stats.csv", "s3://markslab-us-east-2/colabfold/output/popeve/")
    s3_cp_file(f"./{dest_name}_aln_stats.csv", "s3://markslab-private/eve/popeve/data/mappings/")
