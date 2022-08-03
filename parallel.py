#!/usr/bin/env python3
import numpy as np
from multiprocessing import Pool
import pandas as pd

topmed = pd.read_parquet(f'data_singletons_spark_2021_05_25.parquet')
hememap_interactions = pd.read_csv(f'data_hememap_hememap_interactions.csv', sep=',')
hememap_interactions['CHROM'] = hememap_interactions.chr.str.slice(start=3,step=1).astype(str)

A_tot = pd.DataFrame(dict(
        A_id=[1]*5 + [2]*5,
        A_value=range(5, 105, 10)
    ))
B_tot = pd.DataFrame(dict(
        B_id=[1]*5,
        B_low=[0, 30, 30, 46, 84],
        B_high=[10, 40, 50, 54, 84]
    ))

def f(x):
    return x*x

def merge(chr):
    A = A_tot[A_tot.A_id == chr]
    B = B_tot[B_tot.B_id == chr]

    a = A.A_value.values
    aid = A.A_id.values
    bh = B.B_high.values
    bl = B.B_low.values
    bid = B.B_id.values

    i, j = np.where((a[:, None] >= bl) & (a[:, None] <= bh))

    output = pd.concat([
        A.loc[i, :].reset_index(drop=True),
        B.loc[j, :].reset_index(drop=True)
    ], axis=1)

    print(output)

if __name__ == '__main__':
    pools = len(A_tot.A_id.unique())
    with Pool(pools) as p:
        print(p.map(merge, A_tot.A_id.unique()))
