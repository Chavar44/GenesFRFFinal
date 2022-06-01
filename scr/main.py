import os
import pickle
import random
import threading
import time
import yaml
import numpy as np
import pandas as pd
import config
from sklearn.tree import BaseDecisionTree
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


def compute_feature_importance(estimator):
    if isinstance(estimator, BaseDecisionTree):
        return estimator.tree_.compute_feature_importances(normalize=False)
    else:
        importances = [e.tree_.compute_feature_importances(normalize=False)
                       for e in estimator.estimators_]
        importances = np.array(importances)
        return sum(importances, axis=0) / len(estimator)


def separate_dataset(data):
    pass


def train_local_rf(local_data, number_genes):
    # calculate RF/trees for each gene
    # Get the indices of the candidate regulators
    input_idx = list(range(number_genes))
    trees = []
    for i in range(number_genes):
        print('Gene %d/%d...' % (i + 1, number_genes))

        output = local_data[:, i]

        # Normalize output data to unit variance
        # TODO: unit variance must be over whole dataset
        output = output / np.std(output)

        # Remove target gene from candidate regulators
        input_idx = input_idx[:]
        if i in input_idx:
            input_idx.remove(i)

        expr_data_input = local_data[:, input_idx]

        tree_estimator = RandomForestRegressor(n_estimators=config.number_trees, max_features='sqrt')

        # Learn ensemble of trees
        tree_estimator.fit(expr_data_input, output)
        trees.append(tree_estimator)
    return trees


def train(data_hospitals, number_genes):
    # train all local models
    local_random_forrests = []
    for data in data_hospitals:
        local_random_forrests.append(train_local_rf(data, number_genes))

    # train global model
    # TODO: check again for logical errors
    VIM = np.zeros((number_genes, number_genes))
    for i in range(number_genes):
        global_rf = None
        for d in local_random_forrests:
            drf = d[i]

            if global_rf is None:
                global_rf = drf
                global_rf.estimators_ = drf.estimators_
                # global_rf.estimators_ = random.sample(drf.estimators_, trees)
                global_rf.n_estimators = drf.n_estimators
            else:
                global_rf.estimators_ += drf.estimators_
                # global_rf.estimators_ += random.sample(drf.estimators_, trees)
                global_rf.n_estimators += drf.n_estimators
        # compute feature importance for one gene
        feature_importances = compute_feature_importance(global_rf)
        vi = np.zeros(number_genes)
        vi[list(range(number_genes))] = feature_importances
        VIM[i, :] = vi
    VIM = np.transpose(VIM)
    return VIM


if __name__ == "main":
    data = import_data(config.data_path)
    number_patients = data.shape[0]
    number_genes = data.shape[1]
    VIM = np.zeros((number_genes, number_genes))


