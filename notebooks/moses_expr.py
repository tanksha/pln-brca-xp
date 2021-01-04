__author__ = 'Abdulrahman Semrie<hsamireh@gmail.com>'
import pandas as pd
import numpy as np
from mozi_cross_val.main.cross_val import CrossValidation
import os
import argparse
from datetime import datetime

def timer(start_time=None):
    if not start_time:
        start_time = datetime.now()
        return start_time

    elif start_time:
        thour, temp_sec = divmod((datetime.now() - start_time).total_seconds(), 3600)
        tmin, tsec = divmod(temp_sec, 60)
        print('\n Time taken: %i hours %i minutes and %s seconds.' % (thour, tmin, round(tsec, 2)))

def digitize_genes(col, bins=15):
    _, edges = np.histogram(col, bins=bins, density=True)
    res = np.digitize(col, edges)
    return res - res.min()


def digitize_genes_by_median(col):
    median = col.median()
    return col > median


res_ar = np.zeros((16, 4), dtype=np.bool_)
for i in range(16):
    ar0 = bin(i).split('0b')[1]
    ar0 = '0' * (4 - len(ar0)) + ar0
    res_ar[i] = [int(x) for x in ar0]



def binary_genes(merged, genes_columns, by_median=True):
    n_patients, n_genes = merged[genes_columns].shape
    result = dict()
    for gene_col in range(len(genes_columns)):
    # for gene_col in range(1000):
        if by_median:
            digitized = digitize_genes_by_median(merged[genes_columns[gene_col]])
            result[genes_columns[gene_col]] = digitized
        else:
            binary_genes = np.zeros((n_patients, 4), dtype=np.bool_)
            digitized = digitize_genes(merged[genes_columns[gene_col]])
            for i, digit in enumerate(digitized):
                binary_genes[i] = res_ar[digit]
            column_names = [genes_columns[gene_col] + '_{0}'.format(x) for x in range(4)]
            for i, col in enumerate(column_names):
                result[col] = binary_genes[:, i]
    return result



def digitize_non_genes_data(col):
    not_nan = np.nan_to_num(col, nan=-1)
    n_unique = len(set(not_nan))
    if 15 < n_unique:
        n_unique = 15
    edge = np.histogram_bin_edges(np.nan_to_num(col, nan=col.min()), bins=n_unique - 1)
    digits = np.digitize(not_nan, edge)
    digits = digits - digits.min()
    return digits


def preprocess_input(outcome_data, gene_expr_data, target):
    #binarize the gene expression
    outcome_df = pd.read_csv(outcome_data)
    outcome_df = outcome_df[["patient_ID", target]].dropna(axis=0, subset=[target])
    gene_expr_df = pd.read_csv(gene_expr_data)
    gene_cols = gene_expr_df.columns.to_list()[1:]
    patient_ids = gene_expr_df.patient_ID.to_list()
    binary_genes_dict = binary_genes(gene_expr_df, gene_cols)
    binary_genes_dict["patient_ID"] = patient_ids
    result_df = pd.DataFrame(data=binary_genes_dict).sort_values(by='patient_ID') * 1
    merged_df = pd.merge(outcome_df, result_df, on="patient_ID")
    merged_df.drop("patient_ID", axis=1, inplace=True)
    return merged_df


def parse_args():
    parser = argparse.ArgumentParser(description="A script to run moses cross-validation on gene expression data")
    parser.add_argument("--outcome", type=str, default='',
                        help="Path to file that contains patient outcome data")
    parser.add_argument("--expression", type=str, default='',
                        help="Path to file that contains patient gene expression data")
    parser.add_argument("--target_col", type=str, default='posOutcome',
                        help="Name of the target column in the outcome file")
    return parser.parse_args()

def start_moses_run():
    start_time = timer(None)
    print("Starting Moses run")
    args = parse_args()
    target_col = args.target_col
    result_df = preprocess_input(args.outcome, args.expression, target_col)
    moses_opts = "--log-file log1.txt.log --hc-fraction-of-nn 0.01 -j12 --balance 1  --result-count 100 --reduct-knob-building-effort=2  --hc-crossover-min-neighbors=500 --fs-focus=all --fs-seed=init -m 300000 --hc-max-nn-evals=100000 -l debug -q 0.05"
    crossval_options = {"folds": 5, "testSize": 0.3, "randomSeed": 42}
    home_dir = os.environ["HOME"]
    wk_dir = os.path.join(home_dir, "moses_run_" + target_col)
    os.makedirs(wk_dir)
    input_file = result_df.to_csv(os.path.join(wk_dir, "gene_expr_"+target_col+".csv"))

    cross_val = CrossValidation(input_file, wk_dir, target_col, moses_opts, crossval_options, "f1_score", 0.5)

    cross_val.run_folds()
    timer(start_time)
    print("Done")

if __name__ == "__main__":
    start_moses_run()