import numpy as np

from granatum_deep.deep_config import SEED

if SEED:
    np.random.seed(SEED)

from keras.layers import Dense
from keras.layers import Dropout
from keras.layers import Input

from keras.models import Sequential
from keras.models import Model

from keras import regularizers

from time import time

from granatum_deep.deep_config import EPOCHS
from granatum_deep.deep_config import LEVEL_DIMS_IN
from granatum_deep.deep_config import LEVEL_DIMS_OUT
from granatum_deep.deep_config import N_COMPONENTS
from granatum_deep.deep_config import LOSS
from granatum_deep.deep_config import OPTIMIZER
from granatum_deep.deep_config import ACT_REG
from granatum_deep.deep_config import W_REG
from granatum_deep.deep_config import DROPOUT
from granatum_deep.deep_config import ACTIVATION
from granatum_deep.deep_config import DATA_SPLIT


def debug():
    """ """
    X_train = np.random.random((500,10)) * 10
    Y_train = np.random.random((500,5))

    simdeep = DeepBase()
    simdeep.fit(X_train, Y_train)


class DeepBase():
    """ """
    def __init__(self,
                 epochs=EPOCHS,
                 level_dims_in=LEVEL_DIMS_IN,
                 level_dims_out=LEVEL_DIMS_OUT,
                 n_components=N_COMPONENTS,
                 loss=LOSS,
                 optimizer=OPTIMIZER,
                 act_reg=ACT_REG,
                 w_reg=W_REG,
                 dropout=DROPOUT,
                 data_split=DATA_SPLIT,
                 activation=ACTIVATION):
        """
        ### DEFAULT PARAMETER ###:
            dataset=None      ExtractData instance (load the dataset),
            level_dims = [500]
            n_components = 100
            dropout = 0.5
            act_reg = 0.0001
            w_reg = 0.001
            data_split = 0.2
            activation = 'tanh'
            epochs = 10
            loss = 'binary_crossentropy'
            optimizer = 'sgd'
        """
        self.input_dim = None
        self.output_dim = None

        self.epochs = epochs
        self.level_dims_in = level_dims_in
        self.level_dims_out = level_dims_out
        self.n_components = n_components
        self.loss = loss
        self.optimizer = optimizer
        self.dropout = dropout
        self.activation = activation
        self.data_split = data_split

        self.W_l1_constant = w_reg
        self.A_l2_constant = act_reg

        self.encoder_array = {}
        self.model_array = {}

    def predict_autoencoder(self, matrix):
        """
        """
        return self.encoder.predict(matrix)

    def predict_output(self, matrix):
        """
        """
        return self.model.predict(matrix)

    def fit_transform(self, matrix_in, matrix_out=None):
        """
        """
        self.fit(matrix_in, matrix_out)
        return self.transform(matrix_in)

    def fit(self, matrix_in, matrix_out=None):
        """
        main class to create the autoencoder
        """
        if matrix_out is None:
            matrix_out = matrix_in

        self.input_dim = matrix_in.shape[1]
        self.output_dim = matrix_out.shape[1]

        self._create_autoencoder()
        self._compile_model()
        self._fit_autoencoder(matrix_in, matrix_out)

    def _create_autoencoder(self):
        """
        Instantiate the  autoencoder architecture
        """
        t = time()

        model = Sequential()
        nb_hidden = 0

        if hasattr(self, 'n_clusters'):
            self.n_components = self.n_clusters

        ############ Add first side of the AE ######
        for dim in self.level_dims_in:
            nb_hidden += 1
            model = self._add_dense_layer(
                model,
                dim,
                name='hidden layer nb:{0}'.format(nb_hidden))

            if self.dropout:
                model.add(Dropout(self.dropout))

        ############ Add middle layer #############
        model = self._add_dense_layer(
            model,
            self.n_components,
            name='new dim')
        nb_hidden += 1

        if self.dropout:
            model.add(Dropout(self.dropout))

        ############ Add second side of the AE ####
        for dim in self.level_dims_out:
            nb_hidden += 1
            model = self._add_dense_layer(
                model,
                dim,
                name='hidden layer nb:{0}'.format(nb_hidden))

            if self.dropout:
                model.add(Dropout(self.dropout))

        ############ Add final layer ##############
        model = self._add_dense_layer(
            model,
            self.output_dim,
            name='final layer')
        ###########################################
        self.model = model

    def _add_dense_layer(self, model, dim, name=None):
        """
        private function to add one layer
        """
        if not model.layers:
            input_dim = self.input_dim
        else:
            input_dim = None

        model.add(Dense(dim,
                        activity_regularizer=regularizers.l2(self.A_l2_constant),
                        kernel_regularizer=regularizers.l1(self.W_l1_constant),
                        name=name,
                        activation=self.activation,
                        input_dim=input_dim))
        return model

    def _compile_model(self):
        """
        define the optimizer and the loss function
        compile the model and ready to fit the data!
        """
        self.model.compile(optimizer=self.optimizer, loss=self.loss)

    def _fit_autoencoder(self, matrix_in, matrix_out):
        """
        fit the autoencoder using the training matrix
        """

        self.model.fit(x=matrix_in,
                       y=matrix_out,
                       verbose=0,
                       epochs=self.epochs,
                       validation_split=self.data_split,
                       shuffle=True)

        self._define_encoder()

    def _define_encoder(self):
        """ """
        inp = Input(shape=(self.input_dim,))
        encoder = self.model.layers[0](inp)

        if self.model.layers[0].name != 'new dim':

            for layer in self.model.layers[1:]:
                encoder = layer(encoder)
                if layer.name == 'new dim':
                    break

        encoder = Model(inp, encoder)
        self.encoder = encoder

    def transform(self, matrix):
        """
        """
        return self.encoder.predict(matrix)

    def predict(self, matrix):
        """
        """
        return np.argmax(self.transform(matrix), axis=1)

    def fit_predict(self, matrix):
        """
        """
        self.fit(matrix)
        return self.predict(matrix)

    def set_params(self, **kwargs):
        """
        """
        for key in kwargs:
            setattr(self, key, kwargs[key])

        return self


if __name__ == "__main__":
    debug()
