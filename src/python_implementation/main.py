import os.path

import numpy as np
from src.python_implementation import config
from sklearn.tree import BaseDecisionTree
from sklearn.ensemble import RandomForestRegressor
from operator import itemgetter
import multiprocessing as mp


def import_data(path, path_tf=None):
    """
    This function loads the Data from a given path into a numpy array

    :param path: the path to the data given as tsv file
    :param path_tf: if given a path to the transcription factors, they will be loaded as well as the gene_names

    :return data: numpy array with data, gene_names(optional), transcription factors(optional)
    """
    if path_tf is not None:
        data = np.loadtxt(path, dtype=str, skiprows=1)
        raw_gene_names = data[:, :1]
        # Transform raw_gene_names into a list readable by the federated random forest approach
        gene_names = []
        for xs in raw_gene_names:
            for x in xs:
                # Only reads the first 15 characters in order to correctly compare to the Regulators.txt file
                gene_names.append(x[0:15])
        data = data[:, 1:].astype(float)
        tf = np.loadtxt(path_tf, dtype=str)
        return data.T, gene_names, tf.tolist()
    else:
        data = np.loadtxt(path, dtype=str, skiprows=1)[:, 1:].astype(float)
        return data.T


def simulate_different_hospitals(data):
    """
    splits the data into a list of sub data to simulate different hospitals. The number of hospitals is given in the
    config.py file, as well if the dataset is split evenly

    :param data: the data to be split, given as numpy array

    :return: the split data as a list of numpy arrays
    """
    num_pat = data.shape[0]
    data_hospitals = []
    if config.split_even:
        num_pat_per_hospital = num_pat / config.number_of_hospitals
        for i in range(config.number_of_hospitals):
            data_hospitals.append(data[int(num_pat_per_hospital * i): int(num_pat_per_hospital * (i + 1))])
    else:
        assert len(config.split_uneven) == config.number_of_hospitals
        for i in range(config.number_of_hospitals):
            end = len(data) * config.split_uneven[i]
            if i == 0:
                start = 0
            else:
                start = len(data) * config.split_uneven[i - 1]
            data_hospitals.append(data[int(start):int(end)])
    return data_hospitals


def calculate_mean_local(data):
    return np.mean(data, axis=0)


def calculate_std_local(data, global_mean):
    return np.mean((data - global_mean) ** 2, axis=0)


def scaling_of_colums(data, num_genes):
    """
    calculates the std in a federated way, for the scaling to unit variance

    :param data: The data to be scaled (a column)

    :return: the std
    """
    mean_local = np.zeros((len(data), num_genes))
    for i in range(len(data)):
        mean_local[i] = calculate_mean_local(data[i])
    mean_global = np.mean(mean_local, axis=0)

    std_local = np.zeros((len(data), num_genes))
    for i in range(len(data)):
        std_local[i] = calculate_std_local(data[i], mean_global)
    std_global = np.mean(std_local, axis=0)
    return np.sqrt(std_global)


def compute_feature_importance(estimator):
    """
    Computes the feature importance of a tree estimator

    :param estimator: the tree estimator as a trained RandomForestRegressor from sckitLearn

    :return: the mean feature importance of a gene over all patients
    """
    if isinstance(estimator, BaseDecisionTree):
        return estimator.tree_.compute_feature_importances(normalize=False)
    else:
        importances = [e.tree_.compute_feature_importances(normalize=False)
                       for e in estimator.estimators_]
        importances = np.array(importances)
        return np.sum(importances, axis=0) / len(estimator)


def train_local_rf(local_data, std_federated, gene_names=None, regulators='all'):
    """Trains the local random forests for an individual hospital

    Parameters
    ----------
    local_data:
        the local data for one hospital

    std_federated: np.array
        standard variances for each gene to calculate the unit variance

    gene_names: 
        list of strings, optional List of length p, where p is the number of columns in expr_data,
        containing the names of the genes. The i-th item of gene_names must correspond to the i-th column of expr_data.

    regulators:
        list of strings, optional List containing the names of the candidate regulators. When a list of
        regulators is provided, the names of all the genes must be provided (in gene_names). When regulators is set to
        'all', any gene can be a candidate regulator. default: 'all'

    Returns
    -------
        The local trained feature importances
    """

    # TODO: check input
    if not isinstance(local_data, np.ndarray):
        raise ValueError(
            'expr_data must be an array in which each row corresponds to a condition/sample and each column '
            'corresponds to a gene')

    number_genes = local_data.shape[1]

    if gene_names is not None:
        if not isinstance(gene_names, (list, tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != number_genes:
            raise ValueError(
                'input argument gene_names must be a list of length p, where p is the number of columns/genes in the '
                'expr_data')

    if regulators != 'all':
        if not isinstance(regulators, (list, tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            s_intersection = set(gene_names).intersection(set(regulators))
            if not s_intersection:
                raise ValueError('the genes must contain at least one candidate regulator')

    # feature importance matrix
    feature_importance_matrix = np.zeros((number_genes, number_genes))

    for i in range(number_genes):
        print('\tGene %d/%d...' % (i + 1, number_genes))

        output = local_data[:, i]

        if std_federated[i] != 0:
            output = output / std_federated[i]

        # calculation of the indexes to be checked
        if regulators == 'all':
            input_idx = list(range(number_genes))
        else:
            input_idx = [i for i, gene in enumerate(gene_names) if gene in regulators]

        # Remove target gene from candidate regulators
        input_idx = input_idx[:]
        if i in input_idx:
            input_idx.remove(i)

        expr_data_input = local_data[:, input_idx]

        tree_estimator = RandomForestRegressor(n_estimators=config.number_trees, max_features='sqrt')

        # Learn ensemble of trees
        tree_estimator.fit(expr_data_input, output)

        # calculate the feature importance for each gene
        feature_importance = compute_feature_importance(tree_estimator)
        vi = np.zeros(number_genes)
        vi[input_idx] = feature_importance
        # put feature importance into matrix
        feature_importance_matrix[i, :] = vi

    feature_importance_matrix = np.transpose(feature_importance_matrix)

    return feature_importance_matrix


def train(data_hospitals, gene_names=None, regulators='all', parallelize_hospitals=1):
    """Trains each local hospital dataset and calculates the global model

    Parameters
    ----------
    data_hospitals:
        local data from hospitals in a list of numpy arrays

    gene_names: list of strings, optional
        List of length p, where p is the number of columns in expr_data, containing the names of the genes. The i-th item of gene_names must correspond to the i-th column of expr_data.
        default: None

    regulators: list of strings, optional
        List containing the names of the candidate regulators. When a list of regulators is provided, the names of all the genes must be provided (in gene_names). When regulators is set to 'all', any gene can be a candidate regulator.
        default: 'all'

    parallelize_hospitals: int, optional
        Number of processes the individual hospitals run on, must be between 1 and number of different hospitals

    Returns
    -------
        the Feature Importance of the global model
    """
    # TODO: check input!
    if not isinstance(parallelize_hospitals, int):
        raise ValueError(
            'parallelize_hospitals must be an int and strictly positive')
    if 1 > parallelize_hospitals > config.number_of_hospitals:
        raise ValueError('parallelize_hospitals must strictly positive and should not be bigger than the number of '
                         'hospitals')

    # train all local models
    local_feature_importances = []

    number_genes = data_hospitals[0].shape[1]

    std_federated = scaling_of_colums(data_hospitals, number_genes)

    if parallelize_hospitals == 1:
        for index, data in enumerate(data_hospitals):
            file_name = "VIM_H" + str(index + 1) + ".npy"
            path = config.data_path_to_VIM_matrices
            if os.path.exists(os.path.join(path, file_name)):
                print('loading file: ' + file_name)
                local_feature_importances.append(np.load(os.path.join(path, file_name)))
            else:
                print("Hospital %d/%d..." % (index + 1, config.number_of_hospitals))
                local_feature_importances.append(train_local_rf(data, std_federated, gene_names, regulators))
                np.save(os.path.join(path, file_name), local_feature_importances[index])
    else:
        print('running jobs on %d threads' % parallelize_hospitals)
        input_data = list()
        path = config.data_path_to_VIM_matrices
        for i in range(len(data_hospitals)):
            input_data.append([data_hospitals[i], std_federated, gene_names, regulators, i, path])

        pool = mp.Pool(parallelize_hospitals)

        all_output = pool.map(wr_train_local_rf, input_data)

        for (i, vim) in all_output:
            local_feature_importances.append(vim)

    # Calculate the weight of the data of each Hospital
    VIM = np.zeros(local_feature_importances[0].shape)
    if config.split_even:
        for i in range(len(local_feature_importances)):
            VIM += local_feature_importances[i]
        VIM /= len(local_feature_importances)
    else:
        # TODO: Adjust to calculate splits without inside knowledge
        for i in range(0, len(local_feature_importances)):
            VIM = VIM + (local_feature_importances[i] * config.split_uneven[i])
    return VIM


def wr_train_local_rf(args):
    file_name = "VIM_H" + str(args[4] + 1) + ".npy"
    vim_res = train_local_rf(args[0], args[1], args[2], args[3])
    np.save(os.path.join(args[5], file_name), vim_res)
    return [args[4], vim_res]


def get_linked_list_federated(VIM, printing):
    maxcount = 'all'
    ngenes = VIM.shape[0]
    input_idx = range(ngenes)
    vInter = [(i, j, score) for (i, j), score in np.ndenumerate(VIM) if i in input_idx and i != j]

    # Rank the list according to the weights of the edges
    vInter_sort = sorted(vInter, key=itemgetter(2), reverse=True)
    nInter = len(vInter_sort)

    # Random permutation of edges with score equal to 0
    flag = 1
    i = 0
    while flag and i < nInter:
        (TF_idx, target_idx, score) = vInter_sort[i]
        if score == 0:
            flag = 0
        else:
            i += 1

    if not flag:
        items_perm = vInter_sort[i:]
        items_perm = np.random.permutation(items_perm)
        vInter_sort[i:] = items_perm

    # Write the ranked list of edges
    nToWrite = nInter
    if isinstance(maxcount, int) and maxcount >= 0 and maxcount < nInter:
        nToWrite = maxcount

    if printing:
        for i in range(nToWrite):
            (TF_idx, target_idx, score) = vInter_sort[i]
            TF_idx = int(TF_idx)
            target_idx = int(target_idx)
            print('G%d\tG%d\t%.6f' % (TF_idx + 1, target_idx + 1, score))
    return vInter_sort
