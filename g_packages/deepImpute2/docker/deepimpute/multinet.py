import os,binascii

import pandas as pd
import numpy as np
from scipy.stats import pearsonr

import keras
from keras import backend as K
from keras.models import Model,model_from_json
from keras.layers import Dense,Dropout,Input
from keras.callbacks import EarlyStopping
import keras.losses

import tensorflow as tf

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

def generate_random_id():
    rd_number = binascii.b2a_hex(os.urandom(3))

    if type(rd_number) is bytes:
        rd_number = rd_number.decode()
    return rd_number

def get_distance_matrix(raw):

    potential_pred = raw.columns[raw.std() > 0]
    
    covariance_matrix = pd.DataFrame(np.abs(np.corrcoef(raw.T.loc[potential_pred])),
                                     index=potential_pred,
                                     columns=potential_pred).fillna(0)
    return covariance_matrix

def wMSE(y_true,y_pred):
    weights = y_true
    return tf.reduce_mean(weights*tf.square(y_true-y_pred))

class MultiNet:

    def __init__(self,
                 learning_rate=1e-4,
                 batch_size=64,
                 max_epochs=500,
                 patience=5,
                 ncores=20,
                 loss="wMSE",
                 output_prefix="/tmp/multinet",
                 sub_outputdim=512,
                 verbose=1,
                 seed=1234,
                 architecture=None
    ):
        self.NN_parameters = {"learning_rate": learning_rate,
                              "batch_size": batch_size,
                              "loss": loss,
                              "architecture": architecture,
                              "max_epochs": max_epochs,
                              "patience": patience
                              }
        self.sub_outputdim = sub_outputdim
        self.outputdir = "/tmp/{}-{}".format(output_prefix,generate_random_id())
        self.ncores = ncores
        self.verbose = verbose
        self.seed = seed
 
    def loadDefaultArchitecture(self):
        self.NN_parameters['architecture'] = [
                {"type": "dense", "neurons": 256, "activation": "relu"},
                {"type": "dropout", "rate": 0.2},
            ]
        
    def save(self,model):
        os.system("mkdir -p {}".format(self.outputdir))
        
        model_json = model.to_json()
                
        with open("{}/model.json".format(self.outputdir), "w") as json_file:
            json_file.write(model_json)
            
        # serialize weights to HDF5
        model.save_weights("{}/model.h5".format(self.outputdir))
        print("Saved model to disk")

    def load(self):
        json_file = open('{}/model.json'.format(self.outputdir), 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        model = model_from_json(loaded_model_json)
        model.load_weights('{}/model.h5'.format(self.outputdir))

        return model
        
    def build(self,inputdims):
        if self.NN_parameters['architecture'] is None:
            self.loadDefaultArchitecture()

        print(self.NN_parameters['architecture'])

        inputs = [ Input(shape=(inputdim,)) for inputdim in inputdims ]
        outputs = inputs

        for layer in self.NN_parameters['architecture']:
            if layer['type'].lower() == 'dense':
                outputs = [ Dense(layer['neurons'],activation=layer['activation'])(output)
                            for output in outputs ]

            elif layer['type'].lower() == 'dropout':
                outputs = [ Dropout(layer['rate'], seed=self.seed)(output)
                            for output in outputs] 
    
            else:
                print("Unknown layer type.")

        outputs = [ Dense(self.sub_outputdim,activation="softplus")(output)
                    for output in outputs]
                
        model = Model(inputs=inputs,outputs=outputs)

        try:
            loss = eval(self.NN_parameters['loss'])
        except:
            loss = getattr(keras.losses, self.NN_parameters['loss'])
    
        model.compile(optimizer=keras.optimizers.Adam(lr=self.NN_parameters['learning_rate']),
                      loss=loss)

        return model

    def fit(self,
            raw,
            cell_subset=1,
            NN_lim=None,
            genes_to_impute=None,
            ntop=5,
            minVMR=0.5,
            mode='random',
    ):
        if self.seed is not None:
            np.random.seed(self.seed)

        if cell_subset != 1:
            if cell_subset < 1:
                raw = raw.sample(frac=cell_subset)
            else:
                raw = raw.sample(cell_subset)

        gene_metric = (raw.var()/(1+raw.mean())).sort_values(ascending=False)

        if genes_to_impute is None:
            genes_to_impute = self.filter_genes(gene_metric, minVMR, NN_lim=NN_lim)
        else:
            n_genes = len(genes_to_impute)
            if n_genes % self.sub_outputdim != 0:
                print("The number of input genes is not a multiple of {}. Filling with other genes.".format(n_genes))
                fill_genes = gene_metric[:(self.sub_outputdim-n_genes)]
                genes_to_impute = np.concatenate((genes_to_impute,fill_genes))

        covariance_matrix = get_distance_matrix(raw)
        
        self.setTargets(raw.reindex(columns=genes_to_impute), mode=mode)
        self.setPredictors(covariance_matrix,ntop=ntop)

        print("Normalization")
        # normalizer = Normalizer.fromName(self.normalization).fit(raw)
        norm_data = np.log1p(raw).astype(np.float32) # normalizer.transform(raw)

        np.random.seed(self.seed)
        tf.set_random_seed(self.seed)
        
        config = tf.ConfigProto(intra_op_parallelism_threads=self.ncores,
                                inter_op_parallelism_threads=self.ncores,
                                use_per_session_threads=True,
                                allow_soft_placement=True, device_count = {'CPU': self.ncores})
        session = tf.Session(config=config)
        K.set_session(session)

        print("Building network")
        model = self.build([ len(genes) for genes in self.predictors ])

        test_cells = np.random.choice(norm_data.index, int(0.05 * norm_data.shape[0]), replace=False)
        train_cells = np.setdiff1d(norm_data.index, test_cells)

        X_train = [ norm_data.loc[train_cells, inputgenes].values for inputgenes in self.predictors ]
        Y_train = [ norm_data.loc[train_cells, targetgenes].values for targetgenes in self.targets ]
        
        X_test = [ norm_data.loc[test_cells, inputgenes].values for inputgenes in self.predictors ]
        Y_test = [ norm_data.loc[test_cells, targetgenes].values for targetgenes in self.targets ]

        print("Fitting with {} cells".format(norm_data.shape[0]))
        result = model.fit(X_train, Y_train,
                           validation_data=(X_test,Y_test),
                           epochs=self.NN_parameters["max_epochs"],
                           batch_size=self.NN_parameters["batch_size"],
                           callbacks=[EarlyStopping(monitor='val_loss',
                                                    patience=self.NN_parameters["patience"])],
                           verbose=self.verbose)
        print("Stopped fitting after {} epochs".format(len(result.history['loss'])))

        self.save(model)
        
        return self

    def predict(self,
                raw,
                imputed_only=False,
                policy="restore"):

        norm_raw = np.log1p(raw)

        inputs = [ norm_raw.loc[:,predictors].values.astype(np.float32)
                   for predictors in self.predictors ]

        model = self.load()

        predicted = pd.DataFrame(np.hstack(model.predict(inputs)),
                                 index=raw.index,
                                 columns=self.targets.flatten())

        predicted = predicted.groupby(by=predicted.columns,axis=1).mean()
        not_predicted = norm_raw.drop(self.targets.flatten(),axis=1)

        imputed = (pd.concat([predicted,not_predicted],axis=1)
                   .loc[raw.index,raw.columns]
                   .values
        )
        # To prevent overflow
        imputed[ (imputed>2*norm_raw.values.max()) | (np.isnan(imputed)) ] = 0
        # imputed = imputed.mask( (imputed>2*norm_raw.values.max()) | imputed.isnull(), norm_raw)
        # Convert back to counts
        imputed = np.expm1(imputed)
        # imputed = normalizer.transform(imputed,rev=True)        
        
        if policy == "restore":
            print("Filling zeros")
            mask = (raw.values>0)
            imputed[mask] = raw.values[mask]
            # imputed = imputed.mask(raw>0,raw)
        elif policy == "max":
            print("Imputing data with 'max' policy")
            mask = (raw.values>imputed.values)
            imputed[mask] = raw.values[mask]
            # imputed = imputed.mask(raw>imputed,raw)

        imputed = pd.DataFrame(imputed, index=raw.index, columns=raw.columns)

        if imputed_only:
            return imputed.loc[:,predicted.columns]
        else:
            return imputed
        
    def filter_genes(self,
                    gene_metric, # assumes gene_metric is sorted
                    threshold,
                    NN_lim=None
    ):
        if not str(NN_lim).isdigit():
            NN_lim = (gene_metric > threshold).sum()

        n_subsets = int(np.ceil(NN_lim / self.sub_outputdim))
        genes_to_impute = gene_metric.index[:n_subsets*self.sub_outputdim]

        rest = (self.sub_outputdim*n_subsets) % len(genes_to_impute)

        if rest > 0:
            fill_genes = np.random.choice(gene_metric.index, rest)
            genes_to_impute = np.concatenate([genes_to_impute,fill_genes])

        print("{} genes selected for imputation".format(len(genes_to_impute)))

        return genes_to_impute

    def setTargets(self,data, mode='random'):
        
        n_subsets = int(data.shape[1]/self.sub_outputdim)

        if mode == 'progressive':
            self.targets = data.columns.values.reshape([n_subsets,self.sub_outputdim])
        else:
            self.targets = np.random.choice(data.columns,
                                            [n_subsets,self.sub_outputdim],
                                            replace=False)
        
    def setPredictors(self,covariance_matrix,ntop=5):
        
        self.predictors = []
        for i,targets in enumerate(self.targets):
            subMatrix = ( covariance_matrix
                          .loc[targets]
                          .drop(np.intersect1d(targets,covariance_matrix.columns),axis=1)
                          )
            sorted_idx = np.argsort(-subMatrix.values,axis=1)
            predictors = subMatrix.columns[sorted_idx[:,:ntop]].values.flatten()

            self.predictors.append(np.unique(predictors))

            print("Net {}: {} predictors, {} targets"
                  .format(i,len(np.unique(predictors)),len(targets)))

    def score(self,data,policy=None):
        Y_hat = self.predict(data,policy=policy)
        Y = data.loc[Y_hat.index,Y_hat.columns]

        return pearsonr(Y_hat.values.reshape(-1), Y.values.reshape(-1))
        



    
