#!/usr/bin/env python
import os
import subprocess
from multiprocessing import pool
from functools import partial
import collections
import pandas as pd
from tqdm.auto import tqdm


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


tier = 2

if tier == 1:
    s3_cp_file("s3://markslab-us-east-2/colabfold/output/popeve/missing_genes_with_priority.tsv", ".")
    df = pd.read_table("missing_genes_with_priority.tsv")
    assert len(df["protein"].unique()) == len(df)
    prio_df = df[df['high_priority']]

    run_names = [
        # "popeve_priority_proteins_uniref30_2103",
        # "popeve_priority_proteins_uniref30_2202",
        # "popeve_priority_proteins_uniref30_2202_c25",
        "popeve_priority_proteins_uniref30_2202_c40",
        # "popeve_priority_proteins_uniref30_2202_c50",
    ]
if tier == 2:
    s3_cp_file("s3://markslab-us-east-2/colabfold/output/popeve/missing_genes_with_priority.tsv", ".")
    df = pd.read_table("missing_genes_with_priority.tsv")
    assert len(df["protein"].unique()) == len(df)
    prio_df = df[~df['high_priority']]
    run_names = [
        "popeve_tier2_proteins_uniref30_2202_c40",
    ]

def make_a2ms(run_name, tup):
    i, j, protein = tup
    # if subprocess.run(["aws", "s3api", "head-object", "--bucket", "markslab-private", "--key", f"eve/popeve/data/MSA/{run_name}/{protein}.a2m"], capture_output=True).returncode != 0:
    # download a3m
    max_runs = 3
    run = 0
    error = Exception()
    while run < max_runs:
        try:
            subprocess.run(["aws", "s3", "cp", f"s3://markslab-us-east-2/colabfold/output/popeve/{run_name}_pt{i}/{j}.a3m", f"/tmp/{run_name}/{protein}.a3m"], check=True)
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
    os.system("awk '/^>/ { print $0 } !/^>/ { gsub(/[.a-z]/, \"\"); print $1 }' " + f"/tmp/{run_name}/{protein}.a3m > /tmp/{run_name}/{protein}.a2m")
    # export to new location
    run = 0
    while run < max_runs:
        try:
            subprocess.run(["aws", "s3", "cp", f"/tmp/{run_name}/{protein}.a2m", f"s3://markslab-private/eve/popeve/data/MSA/{run_name}/{protein}.a2m"], check=True)
        except subprocess.CalledProcessError as e:
            error = e
            continue
        else:
            break
        finally:
            run += 1
    else:
        raise error
    os.remove(f"/tmp/{run_name}/{protein}.a3m")
    os.remove(f"/tmp/{run_name}/{protein}.a2m")

for run_name in run_names:
    os.makedirs(f"/tmp/{run_name}", exist_ok=True)

    def get_run_files():
        S = 200
        N = int(len(prio_df)/S)
        frames = [prio_df.iloc[i*S:(i+1)*S].copy() for i in range(N+1)]
        for i, frame in enumerate(frames):
            # if i not in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
            #     continue
            for j, row in enumerate(frame.itertuples()):
                yield i, j, row.protein

    p = pool.Pool(len(os.sched_getaffinity(0)) * 2)
    i = tqdm(
        p.imap(
            partial(make_a2ms, run_name),
            get_run_files(),
        ),
        total=len(prio_df)
    )
    collections.deque(i, maxlen=0)

    mapping_df = prio_df[["protein"]].rename(columns={"protein": "protein_name"})
    mapping_df["msa_location"] = mapping_df["protein_name"] + ".a2m"
    mapping_df["theta"] = 0.2
    mapping_df.to_csv(f"./{run_name}_mapping.csv", index=False)
    if not s3_path_exists("markslab-private", f"eve/popeve/data/mappings/{run_name}_mapping.csv"):
        subprocess.run(["aws", "s3", "cp", f"./{run_name}_mapping.csv", "s3://markslab-us-east-2/colabfold/output/popeve/"])
        subprocess.run(["aws", "s3", "cp", f"./{run_name}_mapping.csv", "s3://markslab-private/eve/popeve/data/mappings/"])
