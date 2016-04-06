#!/usr/bin/python
from __future__ import division
import os
#from numpy import *
from scipy import stats
import numpy as np
import operator
import time
import math
import random
import pickle
import copy
from VNF_DPT_source import VNF_request, x, VNF_SPE, CPU_TOTAL, MEM_TOTAL , h_band, brc_base
import matplotlib.pyplot as plt

#flag, indict amp of cpu,mem,band
CPU_flag = 0
Mem_flag = 1
Band_flag = 2
#request num
REQUEST_NUM = 200
#peak num
#Peak_num = 8
#max length of request sequence
MAX_LEN_RS = 20

"""
basically, we mulplex several normal distribution to simulate the time-varying VNF resource request
peak_num is the num of normal distribution that is used.
Every request is a list, each list is consists of seveal VNFs.
and all the request is setted down to a dict. the key of the dicct is the sequence num of the request.
"""

def get_amp(CMB_flag, amp_flag):
    if CMB_flag == CPU_flag:
        if amp_flag > 0.5:
            return random.uniform(2,3)
        else:
            return random.uniform(0.2,0.3)
    if CMB_flag == Mem_flag:
        if amp_flag > 0.5:
            return random.uniform(2,3)
        else:
            return random.uniform(0.2,0.3)
    if CMB_flag == Band_flag:
        if amp_flag > 0.5:
            return random.uniform(2,3)
        else:
            return random.uniform(0.2,0.3)


def get_mu(mu_list):
    mu = random.randint(0,24)
    if len(mu_list) == 0:
        return mu
    else:
        for i_mu in mu_list:
            if i_mu == mu:
                mu = get_mu(mu_list)
        return mu
    

def get_sigma():
    return random.uniform(0.35,0.4)


def get_request_seq():
    request_lib = {}
    VNF_list = []
    for VNF_SPE_i in range(VNF_SPE):
        VNF_list.append(VNF_SPE_i+1)
    print "VNF_list:",VNF_list
    Peak_list = []
    for Peak_i in range(24):
        Peak_list.append(Peak_i)                 #we assume that the peak value arrived at each intiger clock of 24 hours.
    print "Peak_list:", Peak_list
    
    #x = np.arange(0,24,0.1)
    for i in range(REQUEST_NUM):
        species_flag = random.uniform(0,1)
        #if species_flag <= 0.4:          #species_flag <= 0.4 Night_flows, species_flag > 0.9 other flows else, Day_flows
        if i <= 0.5*REQUEST_NUM:
            request_len = random.randint(1,MAX_LEN_RS)
            VNF_slice = random.sample(VNF_list, request_len)
            print "VNF_slice", VNF_slice
            peak_slice = [0,1,2,3,4,5,22,23,24]
            amp_flag = random.uniform(0,1)
            
        #elif species_flag > 0.9:
        elif i > 0.5*REQUEST_NUM and i<= 1.0*REQUEST_NUM:
            request_len = random.randint(1,MAX_LEN_RS)
            VNF_slice = random.sample(VNF_list, request_len)
            print "VNF_slice", VNF_slice
            peak_slice = [8,9,10,11,14,15,16,17,20]
            amp_flag = random.uniform(0,1)
        else:
            request_len = random.randint(1,MAX_LEN_RS)              #generate the length of request randomly, length is the num of VNFs in a request list.
            VNF_slice = random.sample(VNF_list, request_len)        #pick n=request_len different VNFs from the VNF_list randomly
            print "VNF_slice", VNF_slice                            #there is VNF_NUM different vnfs in the simulation 
            peak_num = random.randint(8,12)                         #generate the peak num of the request.
            peak_slice = random.sample(Peak_list, peak_num)
            print "peak_slice:", peak_slice
            amp_flag = random.uniform(0,1)
        #mu = get_mu(mu_list)
        request_list = []
        for j in range(len(VNF_slice)):
            if j == 0:
                VNF_id = VNF_slice[j]
                req_temp = VNF_request()
                for k in range(2):
                    if k == 0:
                        y_cpu_dict = {}
                        y_mem_dict = {}
                        y_band_dict = {}
                        y_cpu = 0
                        y_mem = 0
                        y_band = 0
                        mu_list = []
                        sigma = get_sigma()
                        for mu in peak_slice:
                            y_cpu_dict[mu] = get_amp(0,amp_flag)*stats.norm.pdf(x,mu,sigma) + np.array([float(0.1) for l in range(len(x))])
                            y_mem_dict[mu] = get_amp(1,amp_flag)*stats.norm.pdf(x,mu,sigma) + np.array([float(0.1) for l in range(len(x))])
                            y_band_dict[mu] = get_amp(2,amp_flag)*stats.norm.pdf(x,mu,sigma) + np.array([float(0.1) for l in range(len(x))])
                            y_cpu += y_cpu_dict[mu]
                            y_mem += y_mem_dict[mu]
                            y_band += y_band_dict[mu]                                 
                            mu_list.append(mu)
                        req_temp.id = VNF_id
                        req_temp.cpu = y_cpu
                        req_temp.mem = y_mem
                        req_temp.out_band = y_band
                    if k == 1:
                        y_band_dict = {}
                        y_band = 0
                        mu_list = []
                        sigma = get_sigma()
                        for mu in peak_slice:
                            y_band_dict[mu] = get_amp(2,amp_flag)*stats.norm.pdf(x,mu,sigma) + np.array([float(0.1) for l in range(len(x))])
                            y_band += y_band_dict[mu]                                 
                            mu_list.append(mu)
                        req_temp.in_band = y_band
                request_list.append(req_temp)
            else:
                VNF_id = VNF_slice[j]
                req_temp = VNF_request()
                y_cpu_dict = {}
                y_mem_dict = {}
                y_band_dict = {}
                y_cpu = 0
                y_mem = 0
                y_band = 0
                mu_list = []
                sigma = get_sigma()
                for mu in peak_slice:
                    y_cpu_dict[mu] = get_amp(0,amp_flag)*stats.norm.pdf(x,mu,sigma) + np.array([float(0.1) for l in range(len(x))])
                    y_mem_dict[mu] = get_amp(1,amp_flag)*stats.norm.pdf(x,mu,sigma) + np.array([float(0.1) for l in range(len(x))])
                    y_band_dict[mu] = get_amp(2,amp_flag)*stats.norm.pdf(x,mu,sigma) + np.array([float(0.1) for l in range(len(x))])
                    y_cpu += y_cpu_dict[mu]
                    y_mem += y_mem_dict[mu]
                    y_band += y_band_dict[mu]                                 
                    mu_list.append(mu)
                req_temp.id = VNF_id
                req_temp.cpu = y_cpu
                req_temp.mem = y_mem
                req_temp.in_band = copy.deepcopy(request_list[j-1].out_band)
                req_temp.out_band = y_band
                request_list.append(req_temp)
                #print "mu_list:", mu_list
        if len(request_list) > 1:
            for l in range(len(request_list)):
                if l == 0:
                    request_list[l].next_vnf_id = request_list[l+1].id
                elif l == (len(request_list) - 1):
                    request_list[l].pre_vnf_id = request_list[l-1].id
                else:
                    request_list[l].next_vnf_id = request_list[l+1].id
                    request_list[l].pre_vnf_id = request_list[l-1].id
        
        request_lib[i] = request_list
        
    return request_lib

def search_list_max(list_temp):
    max_value = 0
    list_len = len(list_temp)
    for i in range(list_len):
        if list_temp[i] > max_value:
            max_value = list_temp[i]
    return max_value
    

def check_request(request_lib):       #ensure that the resource of each request to be used less than the host resource,
                            #if the resource of the request more than the host resource, then replace it with the one choosen randomly from the good requests.
    bad_requests = []
    good_requests = []
    for req_id, req_list in request_lib.iteritems():
        req_cpu = 0
        req_mem = 0
        req_band = 0
        for vnf in req_list:
            req_cpu += vnf.cpu
            req_mem += vnf.mem
        req_band = req_list[0].in_band
        cpu_max = search_list_max(req_cpu)
        mem_max = search_list_max(req_mem)
        band_max = search_list_max(req_band)
        if cpu_max > (CPU_TOTAL - len(req_list)*brc_base) or mem_max > (MEM_TOTAL - len(req_list)*brc_base) or band_max > h_band:
            bad_requests.append(req_id)
        else:
            good_requests.append(req_id)
    print "bad_requests",bad_requests, "num:", len(bad_requests)
    for bad_id in bad_requests:
        new_id = random.choice(good_requests)
        request_lib[bad_id] = copy.deepcopy(request_lib[new_id])
    return request_lib

#data = get_request_seq()
#with open('/home/infonet/VNF_deploy_time_varying/data_store/data_new.pickle', 'w') as f:
    #pickle.dump(data,f)
#with open('/home/infonet/VNF_deploy_time_varying/data_store/data_new200.pickle', 'r') as f:
    #data_back = pickle.load(f)
#with open('G:/ldf/simulation/data_store/data_new200.pickle', 'r') as f:
        #data_back = pickle.load(f)

if __name__ == '__main__':
    for i in range(10):
        data_ori = get_request_seq()
        data = check_request(data_ori)
        #with open('/home/infonet/VNF_deploy_time_varying/data_store/data_new500.pickle', 'w') as f:
            #pickle.dump(data,f)
        #with open('/home/infonet/VNF_deploy_time_varying/data_store/data_new500.pickle', 'r') as f:
            #data_back = pickle.load(f)
        
        with open('G:/ldf/simulation/data/'+'data_new200_'+str(i)+'.pickle', 'w') as f:
            pickle.dump(data,f)
        with open('G:/ldf/simulation/data/'+'data_new200_'+str(i)+'.pickle', 'r') as f:
            data_back = pickle.load(f)
    #print "data_back:",data_back
    #x = np.arange(0,24,0.1)
    #data_back_len = len(data_back)
    """
    data_back_keys_list = data_back.keys()
    keys_list_len = len(data_back_keys_list)
    
    
    #for data_i, data_value in data_back.iteritems():
    for i in range(10):
        data_i = data_back_keys_list[i]
        data_value = data_back[data_i]
        print "Request seq_num:",data_i
        plt.figure(data_i)
        VNF_index = 1
        VNF_len = len(data_value)
        for VNF_value in data_value:
            print "data index:", VNF_value.id, "pre_vnf_id", VNF_value.pre_vnf_id, "next_vnf_id", VNF_value.next_vnf_id
            plt.subplot(VNF_len,1,VNF_index)
            plot1, = plt.plot(x,VNF_value.cpu,'cx--')
            plot2, = plt.plot(x,VNF_value.mem,'bo--')
            plot3, = plt.plot(x,VNF_value.in_band,'kp-.')
            plot4, = plt.plot(x,VNF_value.out_band,'rp-.')
            plt.ylabel(VNF_value.id)
            #plt.legend([plot1,plot2,plot3],('cpu','memory','band'),'best',numpoints=1)
            VNF_index += 1
        plt.xlabel('Time')
        plt.legend([plot1,plot2,plot3,plot4],('cpu','memory','in_band','out_band'))#,'best',numpoints=1
    plt.show()
            
            #print "cpu data:", VNF_value.cpu
            #print "mem data:", VNF_value.mem
            #print "band data:", VNF_value.band
    """
    #print "request data:",data



















