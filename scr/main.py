import os
import pickle
import random
import threading
import time
import yaml

import pandas as pd

from distutils import dir_util

from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import train_test_split


def import_data(path):
    """
    import data as Dataframe
    :param path: str
    :return data: pd.Dataframe
    """
    data = pd.read_csv(path)
    return data


def split_data(data):
    test_data = 0
    train_data = 0
    validation_data = 0
    return test_data, train_data, validation_data


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
    pass
