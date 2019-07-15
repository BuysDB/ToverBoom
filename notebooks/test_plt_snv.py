
import pytest
import plot_snv
import pandas as pd
import numpy as np

def test_visual():
    raw = pd.DataFrame.from_dict({'col1': [0, np.nan, 0, np.nan], 'col2': [1, np.nan, np.nan, 1]})
    imputed = pd.DataFrame.from_dict({'col1': [0.45, 0.45, 0.45, 0.45], 'col2': [0.55, 0.55, np.nan, 0.55]})
    combined = pd.DataFrame.from_dict(
        {'col1': {0: 0.0, 1: 0.45, 2: 0.0, 3: 0.45}, 'col2': {0: 1.0, 1: 0.55, 2: np.nan, 3: 1.0}})

    assert plot_snv.gen_imputed_matrix_for_visualization(raw, imputed).equals(combined)

def test_visual2():
    raw = pd.DataFrame.from_dict({'col1': [0, np.nan, 0, np.nan], 'col2': [1, np.nan, np.nan, 1]})
    imputed = pd.DataFrame.from_dict({'col1': [0.45, 0.45, 0.45, 0.45], 'col2': [0.55, 0.55, np.nan, 0.55]})
    combined = pd.DataFrame.from_dict(
        {'col1': {0: 0.0, 1: 0.45, 2: 0.0, 3: 0.45}, 'col2': {0: 1.0, 1: 0.55, 2: np.nan, 3: 1.0}})

    assert plot_snv.gen_imputed_matrix_for_visualization(raw, raw).equals(raw)

def test_construct_df_per_snv_realcase():
    before = pd.read_pickle("./testdata/cellData_before.pickle")
    after = pd.read_pickle("./testdata/cellData_after.pickle")
    assert after.equals(plot_snv.construct_df_per_snv(before, ('chr1', 16352606)))

def test_construct_df_per_snv_fake1():
    before = pd.read_pickle("./testdata/cellData_fake_snv1_5.pickle")
    after = pd.read_pickle("./testdata/cellData_snv1.pickle")
    assert after.equals(plot_snv.construct_df_per_snv(before,'snv1'))
