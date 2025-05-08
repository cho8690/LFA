from colorCheckerConstants import *
import numpy as np

_author_ = "A.J. Aranyosi"
_project_ = "colorChecker"
_file_ = "colorCheckerAffine.py"

"""
(c) 2017-2024 Epicore Biosystems Inc.
Author: A.J. Aranyosi <aja@epicorebiosystems.com>
"""


class ColorCheckerAffine(object):
    """
    TODO: insert class doc
    """

    def __init__(self, parent):
        """
        Constructor for colorCheckerAffine
        """
        self.parent = parent
        self.a_matrix = None
        self.res = None
        self.per_channel_fits = None

    @staticmethod
    def pad(x):
        return np.hstack([x, np.ones((x.shape[0], 1))])

    def get_affine(self, data=None, reference=None):
        if data is None:
            if self.parent.medians is None or not hasattr(self.parent.medians, "__iter__"):
                return False
            data = self.parent.medians
        if reference is None:
            reference = CC_REF_COLORS
        data = np.asarray(data)
        reference = np.asarray(reference)
        try:
            self.a_matrix, self.res, rank, s = np.linalg.lstsq(self.pad(np.asarray(data)), self.pad(np.asarray(reference)),
                                                               rcond=None)
        except TypeError:  # dealing with changes to numpy
            self.a_matrix, self.res, rank, s = np.linalg.lstsq(self.pad(np.asarray(data)), self.pad(np.asarray(reference)))
        return self.a_matrix, self.res

    def apply_affine(self, data, a_matrix=None):
        if a_matrix is None:
            a_matrix = self.a_matrix
        data = np.asarray(data).reshape(-1, a_matrix.shape[1] - 1)
        return np.dot(np.hstack([data, np.ones((data.shape[0], 1))]), a_matrix)[:, :-1]


if __name__ == "__main__":
    colorCheckerAffine = ColorCheckerAffine(None)
