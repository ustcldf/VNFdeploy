#!/usr/bin/python
from __future__ import division
import os
#from numpy import *
from scipy import *
import numpy as np
import operator
import time
import math
from collections import defaultdict
from VNF_DPT_source import VNF_switch, VNF_server, VNF_link, x, root_band, core_band, a_g_band, h_band

"""
Establish the fattree topo
"""

class fat_tree(object):
    def __init__(self, ):
        self.K = 0                       #K determin the size of fat tree
        self.host_num = self.K*self.K*self.K/4
        self.c_s_num = self.K
        self.agg_s_num = self.K*self.K/2
        self.edge_s_num = self.K*self.K/2
        self.root_switch = None
        self.core_id_list = []
        self.agg_id_list = []
        self.edge_id_list = []
        self.core_switch_dict = {}
        self.agg_switch_dict = {}
        self.edge_switch_dict = {}
        self.switch = {}
        self.host_dict = {}
        self.pod_dict = {}
        self.core_link = defaultdict(lambda:defaultdict(lambda:None))
        self.edge_link = defaultdict(lambda:defaultdict(lambda:None))
    
    def pick_even(self, num_list):
        even_list = []
        for num in num_list:
            if num % 2 == 0:
                even_list.append(num)
        even_list.sort()
        return even_list
    
    def pick_odd(self, num_list):
        odd_list = []
        for num in num_list:
            if (num+1) % 2 == 0:
                odd_list.append(num)
        odd_list.sort()
        return odd_list
    
        
    
    
    def instant_topo(self, ):
        #core switch
        self.host_num = self.K*self.K*self.K/4
        self.c_s_num = self.K
        self.agg_s_num = self.K*self.K/2
        self.edge_s_num = self.K*self.K/2
        self.root_switch  = VNF_switch()
        self.root_switch.id = -1
        self.switch[-1] = self.root_switch
        
        core_switch_id = self.K
        for i in range(core_switch_id):
            core_switch = VNF_switch()
            core_switch.id = i
            self.switch[i] = core_switch
            self.core_switch_dict[i] = core_switch
            self.core_id_list.append(i)
        
        end_id = core_switch_id
        for pod_i in range(self.K):
            pod_dict_temp = {}
            agg_list_temp = []
            edge_list_temp = []
            start_id = end_id
            end_id = start_id + self.K
            #agg switch
            for ae_i in range(start_id, int(end_id-self.K/2)):
                agg_switch = VNF_switch()
                agg_switch.id = ae_i
                self.switch[ae_i] = agg_switch
                self.agg_switch_dict[ae_i] = agg_switch
                agg_list_temp.append(ae_i)
                self.agg_id_list.append(ae_i)
            pod_dict_temp['agg'] = agg_list_temp
            #print "self.agg_switch_dict", self.agg_switch_dict
            #edge switch
            for ae_j in range(int(end_id-self.K/2), end_id):
                edge_switch = VNF_switch()
                edge_switch.id = ae_j
                self.switch[ae_j] = edge_switch
                self.edge_switch_dict[ae_j] = edge_switch
                edge_list_temp.append(ae_j)
                self.edge_id_list.append(ae_j)
            pod_dict_temp['edge'] = edge_list_temp
            #store a pod
            self.pod_dict[pod_i] = pod_dict_temp
             
        
        #host
        host_id = int(self.K*self.K*self.K/4)
        for m in range(host_id):
            host = VNF_server()
            host.id = m
            host.band = h_band
            self.host_dict[m] = host
        
        #establish links
        #root links
        for root_i in range(self.K):
            self.core_link[-1][root_i] = np.array([root_band for i in range(len(x))])
            self.core_link[root_i][-1] = np.array([root_band for i in range(len(x))])
        
        
        #core-->agg links
        for c_i in range(int(self.K/2)):
            for a_i in self.pick_even(self.agg_id_list):
                self.core_link[c_i][a_i] = np.array([core_band for i in range(len(x))])
                self.core_link[a_i][c_i] = np.array([core_band for i in range(len(x))])
            
        
        for c_j in range(int(self.K/2),self.K):
            for a_j in self.pick_odd(self.agg_id_list):
                self.core_link[c_j][a_j] = np.array([core_band for i in range(len(x))])
                self.core_link[a_j][c_j] = np.array([core_band for i in range(len(x))])
            
        
        
        for p_i in range(self.K):
            pod_temp = self.pod_dict[p_i]
            for a_g in pod_temp['agg']:
                for e_g in pod_temp['edge']:
                    self.core_link[a_g][e_g] = np.array([a_g_band for i in range(len(x))])
                    self.core_link[e_g][a_g] = np.array([a_g_band for i in range(len(x))])
        end_id_2 = 0
        for p_j in range(self.K):
            pod_temp_2 = self.pod_dict[p_j]
            for e_s in pod_temp_2['edge']:
                start_id_2 = end_id_2
                end_id_2 = start_id_2 + self.K/2
                for h_i in range(int(start_id_2), int(end_id_2)):
                    self.edge_link[e_s][h_i] = np.array([h_band for i in range(len(x))])
        
          
        
    
    
