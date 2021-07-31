# -*- coding: utf-8 -*-
"""
Last modification: 04-03-2020
@author: matthieu.stephant

Tools for python scripts
"""

import numpy as np


def toListInt(x):
    """
    convert np.array or pd.Series into list of int for sending into smart contract
    """
    x = np.array(x*10**6) # to array, and multiplication
    x = np.around(x) # round the numbers
    x = x.astype(int) # from float to int
    x = x.tolist()
    
    return x