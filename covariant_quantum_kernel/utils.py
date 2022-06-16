import pandas as pd
import numpy as np


class QKTCallback:
    """Callback wrapper class."""
    def __init__(self) -> None:
        self._data = [[] for i in range(5)]

    def callback(self, x0, x1=None, x2=None, x3=None, x4=None):
        """
        x[0]: number of function evaluations
        x[1]: the parameters
        x[2]: the function value
        x[3]: the stepsize
        x[4]: whether the step was accepted
        """
        self._data[0].append(x0)
        self._data[1].append(x1)
        self._data[2].append(x2)
        self._data[3].append(x3)
        self._data[4].append(x4)

    def get_callback_data(self):
        return self._data

    def clear_callback_data(self):
        self._data = [[] for i in range(5)]


def train_test_split(data_filepath, num_train=10, num_test=10):
    """Split a dataset into training and test sets."""
    df = pd.read_csv(data_filepath, sep=',', header=None)
    data = df.values

    train = data[:2*num_train, :]
    test = data[2*num_train:2*(num_train+num_test), :]

    ind = np.argsort(train[:,-1])
    X_train = train[ind][:,:-1]
    y_train = train[ind][:,-1]

    ind = np.argsort(test[:,-1])
    X_test = test[ind][:,:-1]
    y_test = test[ind][:,-1]

    return X_train, y_train, X_test, y_test