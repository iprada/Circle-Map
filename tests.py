from __future__ import division

import numpy as np

my_dict = {'A': 0.192, 'T': 0.214, 'C': 0.316, 'G': 0.278}




print(np.sum(value for key, value in my_dict.items() if key != 'A'))