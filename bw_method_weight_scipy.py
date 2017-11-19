#!/usr/bin/env python3
#==========================================================================================
#==========================================================================================
from scipy.optimize import linprog
import numpy as np
from  collections import OrderedDict

def calc_weight(compared2best, compared2worst):
    cb = OrderedDict();
    cw = OrderedDict();
    allkeys = sorted(compared2best.keys());
    for kk in allkeys:
        cb[kk] = compared2best[kk];
        cw[kk] = compared2worst[kk];
    colSize = np.size(allkeys);
    rowSize = 4*colSize-5;
    mat = np.zeros((rowSize-1, colSize+1), dtype=np.double);
    bloc = 0; bkey='';
    wloc = 0; wkey='';
    #print(allkeys)
    # get the best criteria location
    for  keyy, vall in cb.items():
        if vall == 1:
            bkey = keyy;
    bloc = allkeys.index(bkey);
    # get the worst criteria location
    for  keyy, vall in cw.items():
        if vall == 1:
            wkey = keyy;
    wloc = allkeys.index(wkey);
    cb_copy = cb.copy();
    '''delete the key corresponding to best critreia.
    The TWO equations corresponding to the best criterion
    do not appear in the system of equations when comparing all
    criteria against the best '''
    cb_copy.pop(bkey, None);
    tmpmat = np.zeros((len(cb_copy.keys()), colSize+1), dtype=np.double);
    tmpmat1 = np.zeros((len(cb_copy.keys()), colSize+1), dtype=np.double);
    for idx in np.arange(len(cb_copy.keys())):
        itmp = allkeys.index(list(cb_copy.keys())[idx]);
        tmpmat[idx, bloc] = 1.0;
        tmpmat[idx, itmp] = -cb_copy[list(cb_copy.keys())[idx]];
        tmpmat[idx, colSize] = -1.0;
    for idx in np.arange(len(cb_copy.keys())):
        itmp = allkeys.index(list(cb_copy.keys())[idx]);
        tmpmat1[idx, bloc] = -1.0;
        tmpmat1[idx, itmp] = cb_copy[list(cb_copy.keys())[idx]];
        tmpmat1[idx, colSize] = -1.0;
    mat[0:2*colSize-2,:] = np.concatenate((tmpmat, tmpmat1), axis=0);
#----------------------------------------------------------------------------------------
    cw_copy = cw.copy();
    '''delete the two keys corresponding to best critreia and worst criteria.
    The FOUR equations corresponding to the best and worst criteria
    do not appear in the system of equations when comparing all
    criteria against the worst '''
    cw_copy.pop(bkey, None);
    cw_copy.pop(wkey, None);
    tmpmat  = np.zeros((len(cw_copy.keys()), colSize+1), dtype=np.double);
    tmpmat1 = np.zeros((len(cw_copy.keys()), colSize+1), dtype=np.double);
    for idx in np.arange(len(cw_copy.keys())):  
        # find the location of the key of cw_copy in the keylist list
        itmp = allkeys.index(list(cw_copy.keys())[idx]);
        tmpmat[idx, itmp] = 1;
        tmpmat[idx, wloc] = -cw_copy[list(cw_copy.keys())[idx]];
        tmpmat[idx, colSize] = -1.0;
    for idx in np.arange(len(cw_copy.keys())):  
        # find the location of the key of cw_copy in the keylist list
        itmp = allkeys.index(list(cw_copy.keys())[idx]);
        tmpmat1[idx, itmp] = -1;
        tmpmat1[idx, wloc] = cw_copy[list(cw_copy.keys())[idx]];
        tmpmat1[idx, colSize] = -1.0;
    mat[2*colSize-2:,:] = np.concatenate((tmpmat, tmpmat1), axis=0);
    Aeq = np.ones((1, colSize+1), dtype=np.double);
    Aeq[0,-1] = 0.;
    beq = np.array([1]); 
    bub = np.zeros((rowSize-1), dtype=np.double);
    cc = np.zeros((colSize+1), dtype=np.double)
    cc[-1] = 1; 
    res = linprog(cc, A_eq=Aeq, b_eq=beq, A_ub=mat, b_ub=bub, \
                  bounds=(0, None), options={"disp": False});
    sol1 = res['x'];
    outp = dict();
    ii = 0;
    for x in allkeys:
        outp[x] = np.asscalar(sol1[ii]);    
        ii = ii+1;
    return((outp, np.asscalar(sol1[-1])))
   
#==========================================================================================
#==========================================================================================

