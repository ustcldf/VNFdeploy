#!/usr/bin/python
from __future__ import division
import os
import re
import os
import shutil
import sys
from scipy import *
import numpy as np
import operator
import time
import math
import pickle
import copy
from collections import defaultdict

from VNF_dp_algorithm_abs import VNF_simulation_abs
from VNF_dp_algorithm_square import VNF_simulation_square
from VNF_dp_algorithm_cosin import VNF_simulation_cosin
from VNF_dp_algorithm_weighted_cosin import VNF_simulation_w_cosin
from VNF_dp_algorithm_pearson import VNF_simulation_pearson
from VNF_dp_algorithm_xFit import VNF_simulation_xFit_randomized, VNF_simulation_xFit_non_randomized
from VNF_dp_algorithm_xFit_based_vnfs import VNF_simulation_xFit_randomized_w, VNF_simulation_xFit_non_randomized_w
from Algorithm_H import Algorithm_H_r, Algorithm_H_nr


if __name__ == '__main__':
    time_start = time.time()
    for i in range(10):
        with open('G:/ldf/simulation/data/'+'data_new200_'+str(i)+'.pickle', 'r') as f:
            data_back = pickle.load(f)
        
        print "ABS"
        VNF_simulation_abs(data_back)
        
        print "square"
        VNF_simulation_square(data_back)
        
        print "cosin"
        VNF_simulation_cosin(data_back)
        
        print "weighted cosin"
        VNF_simulation_w_cosin(data_back)
        
        print "pearson"
        VNF_simulation_pearson(data_back)
        
        print "xFit_randonm"
        VNF_simulation_xFit_randomized(data_back)
        
        print "xFit_non_randonm"
        VNF_simulation_xFit_non_randomized(data_back)
        
        print "xFit_randonm_w"
        VNF_simulation_xFit_randomized_w(data_back)
        
        print "xFit_non_randonm_w"
        VNF_simulation_xFit_non_randomized_w(data_back)
        
        print "Algorithm H r"
        Algorithm_H_r(data_back)
        
        print "Algorithm H nr"
        Algorithm_H_nr(data_back)
        
        store_path = 'G:/ldf/simulation/data_store/'      #file path, can be changed as needed.
        os.makedirs(store_path+'data_new200_'+str(i))    #file name, can be changed as needed too.
        for subj in os.listdir(store_path):
            #print subj
            #lists = subj.split('.')
            #file_ext = lists[-1]
            #ext_set = ['pickle','json']
            #if os.path.isdir(store_path+subjsubj):
            #if file_ext in ext_set:
            if os.path.isfile(store_path+subj):
                #print "the subj is a file"
                shutil.move(store_path+subj,store_path+'data_new200_'+str(i))
    time_end = time.time()
    time_used = time_end - time_start
    time_second = time_used % 60
    time_ms = (time_used-time_second) / 60
    time_mininute = time_ms % 60
    time_hs = (time_ms - time_mininute) / 60
    time_hour = time_hs % 24
    time_day = (time_hs - time_hour) / 24
    print "used time (seconds)", time_used
    print time_day, "Days",time_hour, "Hours", time_mininute, "Mininutes", time_second, "Seconds"
    print "over"



