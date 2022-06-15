import numpy as np
from src.python_implementation import config
from sklearn.tree import BaseDecisionTree
from sklearn.ensemble import RandomForestRegressor
from operator import itemgetter



def import_data(path):
    """
    This function loads the Data from a given path into a numpy array

    :param path: the path to the data given as tsv file

    :return data: numpy array with data
    """
    data = np.loadtxt(path, dtype=str, skiprows=1)[:, 1:].astype(float)
    return data


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
            data_hospitals.append(data[int(num_pat_per_hospital*i): int(num_pat_per_hospital*(i+1))])
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


def scaling_of_colums(data):
    """
    Scales the column of a column (gene) to unit variance over all patients

    :param data: The data to be scaled (a column)

    :return: the scaled data
    """
    # TODO: ImplementMe
    # scale the columns with unit variance
    pass


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


def train_local_rf(local_data, number_genes):
    """
    Trains the local random forrests for a individual hospital

    :param local_data: the local data for one hospital
    :param number_genes: the number of genes in the dataset

    :return: Either!!!! The local tree or the the local trained feature importances
    """
    # calculate RF/trees for each gene
    # Get the indices of the candidate regulators

    trees = []
    for i in range(number_genes):
        print('\tGene %d/%d...' % (i + 1, number_genes))

        input_idx = list(range(number_genes))
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
    """
    Trains each local hospital dataset and calculates the global model

    :param data_hospitals: local data from hospitals in a list of numpy arrays
    :param number_genes: the number of genes in the dataset

    :return: the Feature Importance of the global model
    """
    # train all local models
    local_random_forrests = []
    for index, data in enumerate(data_hospitals):
        print("Hospital %d/%d..." % (index + 1, config.number_of_hospitals))
        local_random_forrests.append(train_local_rf(data, number_genes))

    # train global model
    # TODO: check again for logical errors
    VIM = np.zeros((number_genes, number_genes))
    for i in range(number_genes):
        # TODO: remove instead have local feature importances
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
        # TODO: move to local calculations
        feature_importances = compute_feature_importance(global_rf)
        # TODO: calculate mean over all local feature importances
        vi = np.zeros(number_genes)
        input_idx = list(range(number_genes))
        if i in input_idx:
            input_idx.remove(i)
        vi[input_idx] = feature_importances
        VIM[i, :] = vi
    VIM = np.transpose(VIM)
    return VIM


def get_linked_list_federated(VIM):
    file_name = None
    gene_names = None
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

    if file_name:

        outfile = open(file_name, 'w')

        if gene_names is not None:
            for i in range(nToWrite):
                (TF_idx, target_idx, score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                outfile.write('%s\t%s\t%.6f\n' % (gene_names[TF_idx], gene_names[target_idx], score))
        else:
            for i in range(nToWrite):
                (TF_idx, target_idx, score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                outfile.write('G%d\tG%d\t%.6f\n' % (TF_idx + 1, target_idx + 1, score))

        outfile.close()

    else:

        if gene_names is not None:
            for i in range(nToWrite):
                (TF_idx, target_idx, score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                print('%s\t%s\t%.6f' % (gene_names[TF_idx], gene_names[target_idx], score))
        else:
            for i in range(nToWrite):
                (TF_idx, target_idx, score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                print('G%d\tG%d\t%.6f' % (TF_idx + 1, target_idx + 1, score))


if __name__ == '__main__':
    data = import_data(config.data_path)
    number_patients = data.shape[0]
    number_genes = data.shape[1]
    hospital_data = simulate_different_hospitals(data)
    vim = train(hospital_data, number_genes)
    print(vim)



