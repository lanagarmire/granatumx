import os
import numpy as np
import pandas as pd
import binascii
import warnings
import tempfile
from math import ceil
from multiprocessing import cpu_count, sharedctypes
from multiprocessing.pool import Pool
from sklearn.metrics import r2_score

from deepimpute.net import Net
from deepimpute.normalizer import Normalizer
from deepimpute.util import get_input_genes,get_target_genes
from deepimpute.util import score_model

def newCoreInitializer(arr_to_populate):
    global sharedArray
    sharedArray = arr_to_populate

def trainNet(in_out, NN_param_i, data_i, labels):
    features, targets = in_out

    net = Net(**NN_param_i)
    net.fit(data_i, targetGenes=targets, predictorGenes=features, labels=labels)

    # retrieve the array
    params = list(NN_param_i.keys()) + ['targetGenes', 'NNid', 'predictorGenes']
    args2return = [(attr, getattr(net, attr)) for attr in params]
    return {k: v if k[0] != '_' else (k[1:], v) for k, v in args2return}

def predictNet(data_i, NN_param_i, labels):
    net = Net(**NN_param_i)
    data_i_ok = pd.DataFrame(np.reshape(data_i, list(map(len, labels))),
                             index=labels[0], columns=labels[1])
    return net.predict(data_i_ok)

def trainOrPredict(args):
    in_out, NN_param_i, labels, mode = args
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        data_i = np.ctypeslib.as_array(sharedArray)
    if mode == "predict":
        return predictNet(data_i, NN_param_i, labels)
    return trainNet(in_out, NN_param_i, data_i, labels)


class MultiNet(object):
    def __init__(self, n_cores=4, predictorLimit=10, preproc='log_or_exp', runDir=os.path.join(tempfile.gettempdir(),'run'), seed=0, **NN_params):
        self._maxcores = n_cores
        self.predictorLimit = predictorLimit
        self.norm = Normalizer.fromName(preproc)
        self.runDir = runDir
        self.seed = seed
        self.NN_params = NN_params
        self.seed = seed
        self.NN_params['seed'] = seed

        if 'dims' not in self.NN_params.keys():
            self.NN_params['dims'] = [20,500]

    @property
    def maxcores(self):
        if self._maxcores == 'all':
            return cpu_count()
        else:
            return self._maxcores

    @maxcores.setter
    def maxcores(self, value):
        self._maxcores = value

    def get_params(self, deep=False):
        return self.__dict__

    def setIDandRundir(self,data):
        # set runID
        runID = binascii.b2a_hex(os.urandom(5))
        if type(runID) is bytes:
            runID = runID.decode()
        self.NN_params['runDir'] = os.path.join(self.runDir, str(runID))

    def getCores(self,NN_genes):
        n_runs = int(ceil(1.*len(NN_genes) / self.NN_params['dims'][1]))
        n_cores = min(self.maxcores, n_runs)
        self.NN_params['n_cores'] = max(1, int(self.maxcores / n_cores))
        return n_runs,n_cores

    def fit(self, data, NN_lim='auto', cell_subset=None):
        np.random.seed(seed=self.seed)

        df = pd.DataFrame(data)

        self.setIDandRundir(df)

        # Change the output dimension if the data has too few genes
        if df.shape[1] < self.NN_params['dims'][1]:
            self.NN_params['dims'][1] = df.shape[1]

        # Choose genes to impute
        genes_sort = df.quantile(.99).sort_values(ascending=False)
        NN_genes = get_target_genes(genes_sort,NN_lim=NN_lim)

        df_to_impute = df[NN_genes]

        n_runs,n_cores = self.getCores(NN_genes)


        # ------------------------# Subnetworks #------------------------#

        predictors = np.intersect1d(genes_sort.index[genes_sort>self.predictorLimit], NN_genes)
        print('Using {} genes as potential predictors'.format(len(predictors)))

        n_choose = int(len(NN_genes)/self.NN_params['dims'][1])

        subGenelists = np.random.choice(NN_genes,
                                        [n_choose, self.NN_params['dims'][1]],
                                        replace=False).tolist()
        if n_choose < n_runs:
            # Special case: for the last run, the output layer will have less nodes
            selectedGenes = np.reshape(subGenelists, -1)
            subGenelists.append(np.setdiff1d(NN_genes, selectedGenes).tolist())

        # ------------------------# Extracting input genes #------------------------#

        corrMatrix = 1 - np.abs(pd.DataFrame(np.corrcoef(df_to_impute.T),
                                             index=NN_genes, columns=NN_genes)[predictors])

        in_out_genes = get_input_genes(df_to_impute,self.NN_params['dims'],distanceMatrix=corrMatrix,
                                       targets=subGenelists,predictorLimit=self.predictorLimit)

        # ------------------------# Subsets for fitting #------------------------#

        n_cells = df_to_impute.shape[0]

        if type(cell_subset) is float or cell_subset == 1:
            n_cells = int(cell_subset * n_cells)

        elif type(cell_subset) is int:
            n_cells = cell_subset

        self.trainCells = df_to_impute.sample(n_cells,replace=False).index

        print('Starting training with {} cells ({:.1%}) on {} threads ({} cores/thread).'.
              format(n_cells, 1.*n_cells/df_to_impute.shape[0], n_cores, self.NN_params['n_cores']))

        # -------------------# Preprocessing (if any) #--------------------#

        df_to_impute = self.norm.fit(df_to_impute).transform(df_to_impute)

        # -------------------# Share matrix between subprocesses #--------------------#

        ''' Create memory chunk and put the matrix in it '''
        idx, cols = self.trainCells, df_to_impute.columns
        trainData = df_to_impute.loc[self.trainCells, :].values

        ''' Parallelize process with shared array '''
        childJobs = [(in_out, self.NN_params, (idx, cols), 'train')
                      for in_out in in_out_genes]

        output_dicts = self.runOnMultipleCores(n_cores, trainData.flatten(), childJobs)

        self.networks = []
        for dictionnary in output_dicts:
            self.networks.append(Net(**dictionnary))

        return self

    def runOnMultipleCores(self, cores, data, childJobs):
        sharedArray = sharedctypes.RawArray('d', data)

        pool = Pool(processes=cores, initializer=newCoreInitializer, initargs=(sharedArray,))
        output_dicts = pool.map(trainOrPredict, childJobs)
        pool.close()
        pool.join()
        return output_dicts


    def predict(self, data, imputed_only=False, restore_pos_values=True):

        df = pd.DataFrame(data)

        ''' Create memory chunk and put the matrix in it '''
        idx, cols = df.index, df.columns
        df_norm = self.norm.fit(df).transform(df).values.flatten()

        ''' Parallelize process with shared array '''
        childJobs = [((12, 15), net.__dict__, (idx, cols), 'predict')
                  for net in self.networks]

        output_dicts = self.runOnMultipleCores(self.maxcores, df_norm, childJobs)

        Y_imputed = pd.concat(output_dicts, axis=1)
        Y_not_imputed = df[[gene for gene in df.columns if gene not in Y_imputed.columns]]
        Y_total = self.norm.transform(pd.concat([Y_imputed, Y_not_imputed], axis=1)[df.columns],
                                      rev=True)
        if restore_pos_values:
            Y_total = Y_total.mask(df>0,df)
        if imputed_only:
            Y_total = Y_total[Y_imputed.columns]

        if type(data) == type(pd.DataFrame()):
            return Y_total
        else:
            return Y_total.values

    def score(self, data, metric=r2_score):
        imputedGenes = list(zip(*[ net.targetGenes for net in self.networks ]))
        return score_model(self,pd.DataFrame(data),metric=r2_score, cols=imputedGenes)
