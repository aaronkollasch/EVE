#!/usr/bin/env python
import os
import subprocess
from multiprocessing import pool
from functools import partial
import collections
import pandas as pd
from tqdm.auto import tqdm
from Bio.SeqIO.FastaIO import SimpleFastaParser


def s3_path_exists(bucket, key, silent=False):
    result = subprocess.run(
        ["aws", "s3api", "head-object", "--bucket", bucket, "--key", key],
        stdout=subprocess.DEVNULL if silent else None,
        stderr=subprocess.DEVNULL if silent else None,
    )
    if result.returncode == 0:
        return True
    elif result.returncode == 255:
        return False
    else:
        raise Exception(result.returncode, result.stderr)

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


tier = 1

if tier == 1:
    project_name = "proteingym"
    df = pd.read_csv("clinical_subs_need_weights_2023_10_17_EVE.csv")
    df.rename(columns={"protein_name": "protein"}, inplace=True)
    assert len(df["protein"].unique()) == len(df)
    prio_df = df
    run_names = [
        "proteingym_colabfold_to_run_uniref30_2202_c50",
    ]
    prio_df["protein_index"] = -1
    with open("colabfold_to_run.fasta") as f:
        for i, (name, seq) in enumerate(SimpleFastaParser(f)):
            name = name.split()[0]
            prio_df.loc[prio_df["protein"] == name, "protein_index"] = i
    prio_df = prio_df[(prio_df["protein_index"] != -1) & (prio_df["msa_location"].str.contains("_colabfold"))]
    print(len(df), len(prio_df))
    print(prio_df)

if tier == 2:
    project_name = "proteingym"
    df = pd.read_csv("clinical_subs_need_weights_2023_10_17_EVE.csv")
    df.rename(columns={"protein_name": "protein"}, inplace=True)
    assert len(df["protein"].unique()) == len(df)
    prio_df = df
    run_names = [
        "proteingym_clinical_subs_uniref30_2202_c50_20231018",
    ]
    prio_df["protein_index"] = -1
    with open("clinical_subs_need_weights_2023_10_18.fasta") as f:
        for i, (name, seq) in enumerate(SimpleFastaParser(f)):
            name = name.split()[0]
            prio_df.loc[prio_df["protein"] == name, "protein_index"] = i
    prio_df = prio_df[(prio_df["protein_index"] != -1) & (prio_df["msa_location"].str.contains("_colabfold"))]
    print(len(df), len(prio_df))
    print(prio_df)

for name in run_names:
    os.makedirs(f"/tmp/{name}/", exist_ok=True)

def make_a2ms(run_name, project_name, tup):
    i, msa_location, protein = tup
    msa_location = msa_location.replace(".a2m", "")
    # download a3m
    max_runs = 3
    run = 0
    error = Exception()
    while run < max_runs:
        try:
            subprocess.run(["aws", "s3", "cp", f"s3://markslab-us-east-2/colabfold/output/{project_name}/{run_name}/{i}.a3m", f"/tmp/{run_name}/{msa_location}.a3m"], check=True)
        except subprocess.CalledProcessError as e:
            error = e
            continue
        else:
            break
        finally:
            run += 1
    else:
        raise error
    # convert to a2m
    os.system("awk '/^>/ { print $0 } !/^>/ { gsub(/[.a-z]/, \"\"); print $1 }' " + f"/tmp/{run_name}/{msa_location}.a3m > /tmp/{run_name}/{msa_location}.a2m")
    # export to new location
    run = 0
    while run < max_runs:
        try:
            subprocess.run(["aws", "s3", "cp", f"/tmp/{run_name}/{msa_location}.a2m", f"s3://markslab-private/eve/{project_name}/data/MSA/{run_name}/{msa_location}.a2m"], check=True)
        except subprocess.CalledProcessError as e:
            error = e
            continue
        else:
            break
        finally:
            run += 1
    else:
        raise error
    os.remove(f"/tmp/{run_name}/{msa_location}.a3m")
    os.remove(f"/tmp/{run_name}/{msa_location}.a2m")

for run_name in run_names:
    os.makedirs(f"/tmp/{run_name}", exist_ok=True)

    def get_run_files():
        for row in prio_df.itertuples():
            yield row.protein_index, row.msa_location, row.protein

    p = pool.Pool(len(os.sched_getaffinity(0)) * 2)
    #p = pool.Pool(1)
    i = tqdm(
        p.imap(
            partial(make_a2ms, run_name, project_name),
            get_run_files(),
        ),
        total=len(prio_df)
    )
    collections.deque(i, maxlen=0)

    # mapping_df = prio_df[["protein"]].rename(columns={"protein": "protein_name"})
    # mapping_df["msa_location"] = mapping_df["protein_name"] + ".a2m"
    # mapping_df["theta"] = 0.2
    # mapping_df.to_csv(f"./{run_name}_mapping.csv", index=False)
    # if not s3_path_exists("markslab-private", f"eve/{project_name}/data/mappings/{run_name}_mapping.csv"):
    #     subprocess.run(["aws", "s3", "cp", f"./{run_name}_mapping.csv", "s3://markslab-us-east-2/colabfold/output/{project_name}/"])
    #     subprocess.run(["aws", "s3", "cp", f"./{run_name}_mapping.csv", "s3://markslab-private/eve/{project_name}/data/mappings/"])
