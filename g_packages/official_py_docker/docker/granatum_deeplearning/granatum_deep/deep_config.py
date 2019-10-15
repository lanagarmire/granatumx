"""
FILE TO PUT GLOBAL VARIABLES AND PARAMETERS
"""

##################### Autoencoder Variable ##############
## Dimensions of the intermediate layers before and after the middle hidden layer
# level on the dims of the hidden layers BEFORE new dim
LEVEL_DIMS_IN = [50]
# Number of nodes in the middle hidden layer
# (i.e. the new dimensions of the transformed data)
N_COMPONENTS = 100
# level on the dims of the hidden layers AFTER new dim
LEVEL_DIMS_OUT = [50]
# Percentage of edges being dropout at each training iteration (None for no dropout)
DROPOUT = 0.5
# L2 Regularization constant on the node activity
ACT_REG = False
# L1 Regularization constant on the weight
W_REG = False
# Fraction of the dataset to be used as test set when building the autoencoder
DATA_SPLIT = None
# activation function
ACTIVATION = 'relu'
# Number of epoch
EPOCHS = 10
# Loss function to minimize
LOSS = 'binary_crossentropy'
# Optimizer (sgd for Stochastic Gradient Descent)
OPTIMIZER = 'adam'
# train/test split
DATA_SPLIT = None
# Random seed
SEED = None
##########################################################
