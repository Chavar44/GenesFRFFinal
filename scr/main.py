import os
import pickle
import random
import threading
import time
import yaml
import numpy as np
import pandas as pd
import config

from distutils import dir_util

from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import train_test_split


def import_data(path):
    """
    import data as Dataframe
    :param path: str
    :return data: pd.Dataframe
    """
    data = np.loadtxt(path, dtype=str, skiprows=1)[:, 1:].astype(float)
    return data


def split_data(data):
    test_data = 0
    train_data = 0
    validation_data = 0
    return test_data, train_data, validation_data


def simulate_different_hospitals(data, num_pat):
    data_hospitals = []
    if config.split_even:
        num_pat_per_hospital = num_pat/config.number_of_hospitals
        for i in range(config.number_of_hospitals):
            data_hospitals.append(data[num_pat_per_hospital*i: num_pat_per_hospital*(i+1)])
    else:
        # TODO: create unbalanced data set
        pass
    return data_hospitals


def scaling_of_colums(data):
    # scale the colums with unit variance
    pass


def feature_importance(data):
    pass


def separate_dataset(data):
    pass


def train_local_rf():
    pass


def train():
    pass


def predict():
    pass


if __name__ == "main":
    data = import_data(config.data_path)
    number_patients = data.shape[0]
    number_genes = data.shape[1]
    VIM = np.zeros((number_genes, number_genes))


