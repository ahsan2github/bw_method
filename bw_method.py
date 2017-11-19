#!/usr/bin/env python3
from bw_method_weight_scipy import calc_weight
import numpy as np

compared2best  = dict({'price':1, 'mp':2, 'fps':7,  'ISO':4, 'video':8, 'lenseco':2, 'conn':5});
compared2worst = dict({'price':9, 'mp':5, 'fps':8,  'ISO':1, 'video':1, 'lenseco':4, 'conn':5});

w, zeta = calc_weight(compared2best, compared2worst)
print("Weights: ", w);
print("sum(Weights): ", np.sum(list(w.values())));
print('Consistency: {0}'.format(zeta));



