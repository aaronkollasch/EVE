#!/usr/bin/env python
import os
import subprocess
from multiprocessing import pool
from functools import partial
import collections
import pandas as pd
from tqdm.auto import tqdm

df = pd.read_table("missing_genes_with_priority.tsv")
assert len(df["protein"].unique()) == len(df)
prio_df = df[df['high_priority']]

# mapping_df = prio_df[["protein"]].rename(columns={"protein": "protein_name"})
# mapping_df["msa_location"] = mapping_df["protein_name"] + ".a2m"
# mapping_df["theta"] = 0.2
# mapping_df.to_csv(f"./{run_name}_mapping.csv", index=False)
# subprocess.run(["aws", "s3", "cp", f"./{run_name}_mapping.csv", "s3://markslab-us-east-2/colabfold/output/popeve/"])
# subprocess.run(["aws", "s3", "cp", f"./{run_name}_mapping.csv", "s3://markslab-private/eve/popeve/data/mappings/"])

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

for run_name in [
#     "popeve_priority_proteins_uniref30_2103",
    # "popeve_priority_proteins_uniref30_2202",
    # "popeve_priority_proteins_uniref30_2202_c25",
    "popeve_priority_proteins_uniref30_2202_c40",
    # "popeve_priority_proteins_uniref30_2202_c50",
]:
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
