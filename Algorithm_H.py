#!/usr/bin/python
from __future__ import division
import os
#from numpy import *
from scipy import *
import numpy as np
import operator
import time
import math
import copy
import json
import pickle
import random
from collections import defaultdict
from VNF_topo_generate import fat_tree
#from VNF_request_generate import REQUEST_NUM, get_request_seq, data_back
from VNF_DPT_source import x, VNF_SPE, brc_base, cpu_brc_per, mem_brc_per, CPU_TOTAL,MEM_TOTAL,h_band,band_alpha

ROOT_ID = -1
INFINITY = 1e10
EPSILON_F = 5.960465e-8

class VNF_route(object):
    def __init__(self, ):
        self.switches = None
        self.adjacency = None #defaultdict(lambda:defaultdict(lambda:None))
        # [sw1][sw2] -> (distance, intermediate)
        self.path_map = defaultdict(lambda:defaultdict(lambda:(None,None)))
    
    def _calc_paths(self, ):
        def dump ():
          for i in sws.iteritems():
            for j in sws.iteritems():
              a = self.path_map[i][j][0]
              #a = adjacency[i][j]
              if a is None: a = "*"
              print a
            # print
      
        sws = self.switches
        self.path_map.clear()
        for k,s in sws.iteritems():
          for j,value in self.adjacency[s.id].iteritems():
            if value is None: continue
            self.path_map[s.id][j] = (1,None)
          self.path_map[s.id][s.id] = (0,None) # distance, intermediate
        #dump()
      
        for k,s_1 in sws.iteritems():
          for i,s_2 in sws.iteritems():
            for j,s_3 in sws.iteritems():
              if self.path_map[s_2.id][s_1.id][0] is not None:
                if self.path_map[s_1.id][s_3.id][0] is not None:
                  # i -> k -> j exists
                  ikj_dist = self.path_map[s_2.id][s_1.id][0]+self.path_map[s_1.id][s_3.id][0]
                  if self.path_map[s_2.id][s_3.id][0] is None or ikj_dist < self.path_map[s_2.id][s_3.id][0]:
                    # i -> k -> j is better than existing
                    self.path_map[s_2.id][s_3.id] = (ikj_dist, s_1.id)
      
        #print "--------------------"
        #dump()
    
    def _get_raw_path (self, src, dst):
        """
        Get a raw path (just a list of nodes to traverse)
        """
        if len(self.path_map) == 0: self._calc_paths()
        if src is dst:
          # We're here!
          return []
        if self.path_map[src][dst][0] is None:
          return None
        intermediate = self.path_map[src][dst][1]
        if intermediate is None:
          # Directly connected
          return []
        return self._get_raw_path(src, intermediate) + [intermediate] + \
               self._get_raw_path(intermediate, dst)
    
    def _check_path (self, p):
        """
        Make sure that a path is actually a string of nodes
      
        returns True if path is valid
        """
        for a,b in zip(p[:-1],p[1:]):
          if self.adjacency[a][b] == None:
            return False
          if self.adjacency[b][a] == None:
            return False
        return True
    
    def _get_path (self, src, dst):
        """
        Gets a cooked path -- a list of (node)
        """
        # Start with a raw path...
        if src == dst:
          path = [src]
        else:
          path = self._get_raw_path(src, dst)
          if path is None: return None
          path = [src] + path + [dst]
        #print "what's src is:",src
        #print "what's dst is:",dst
        #print "what's in path is:",'\n',path
      
        assert self._check_path(path), "Illegal path!"
      
        return path
    
    def _get_switch(self,host_id,e_link):
        for id_1,iterm_1 in e_link.iteritems():
            for id_2,iterm_2 in iterm_1.iteritems():
                if id_2 == host_id:
                    return id_1


class VNF_algorithm_h(object):
    
    def verify_node_resource_left(self, path, link_map, request, host):
        #sw_id = route._get_switch(host.id)
        #path = route._get_path(ROOT_ID, sw_id)
        #modify host resource left
        req_len = len(request)               #here, request is a list
        req_cpu = 0
        req_mem = 0
        req_band = 0
        req_cpu += request[1].cpu
        req_mem += request[1].mem
        
        #modify the vnfs in each host
        vnf_in_req = []
        vnf_in_req.append(request[1].id)
        #print "vnf num in req", len(vnf_in_req), vnf_in_req
        vnf_add_new = []
        for vnf_i in vnf_in_req:
            if vnf_i in host.vnf_container:
                continue
            else:
                host.vnf_container.append(vnf_i)
                vnf_add_new.append(vnf_i)
                host.cpu_res -= cpu_brc_per
                host.mem_res -= mem_brc_per
        
        host.verify_cpu(req_cpu)
        host.verify_mem(req_mem)
    
    def check_value(self, list_ori, list_cmp):             #our data is a sequence, check value is to calc the break times of resource
        list_len = len(list_cmp)                #list_ori is the left resource sequence, list_cmp is the resource sequence that is to be setted to responding host
        break_t = 0
        for i in range(list_len):
            if list_ori[i] - list_cmp[i] < 0 :
                break_t += 1
        return break_t
        
    def check_positive(self, list_temp):
        list_len = len(list_temp)
        positive_num = 0
        for i in range(list_len):
            if list_temp[i] < 0:
                positive_num += 1
        return positive_num
        
    def calc_req_vnf(self, req_list):                #count the vnf species in req
        vnf_list = []
        req_len = len(req_list)
        for i in range(req_len):
            vnf_list.append(req_list[i].id)
        
        return vnf_list

    def check_end(self, path, link_map, request, host):          #judge if the request can be setted to host
        #host_backup = copy.deepcopy(host)
        req_len = len(request)             #request is a list
        req_cpu = 0
        req_mem = 0
        req_band = 0
        req_cpu = request[1].cpu
        req_mem = request[1].mem
        req_band = request[1].in_band + request[1].out_band                           #request band is determined by the first VNF
        vnf_species = request[1].id
        #print "vnf num in req", len(vnf_in_req), vnf_in_req
        vnf_add_new = []
        if vnf_species in host.vnf_container:
            pass
        else:
            host.vnf_container.append(vnf_species)
            vnf_add_new.append(vnf_species)
            host.cpu_res -= cpu_brc_per
            host.mem_res -= mem_brc_per
        host_vnf_container = copy.deepcopy(host.vnf_container)
        sw_band_bt_total = 0
        
        for a,b in zip(path[:-1],path[1:]):
            sw_band_bt_total  += self.check_value(link_map[a][b], req_band)
            #band_b_t_list.append(band_b_t)
        h_cpu_bt_total = self.check_value(host.cpu_res, req_cpu)                       #break times
        h_mem_bt_total = self.check_value(host.mem_res, req_mem)
        h_band_bt_total = self.check_value(host.band_res, req_band)
        
        
        #print "h_cpu_bt_total, h_mem_bt_total, h_band_bt_total",h_cpu_bt_total, h_mem_bt_total, h_band_bt_total
        if sw_band_bt_total == 0 and h_cpu_bt_total == 0 and h_mem_bt_total == 0 and h_band_bt_total == 0:                     #to be changed, a threshold can be defined
            #print "resource verified!"
            
            if len(vnf_add_new) != 0:
                for vnf_j in vnf_add_new:         #just check, no matter the host could hold the request or not, it'll should be restore back
                    host.vnf_container.remove(vnf_j)
                    host.cpu_res += cpu_brc_per
                    host.mem_res += mem_brc_per
                #host_vnf_container_2 = copy.deepcopy(host.vnf_container)
            
            return 1
        else:
            #host = copy.deepcopy(host_backup)
            if len(vnf_add_new) != 0:
                for vnf_k in vnf_add_new:         #if can not hold, restore the vnf_container
                    host.vnf_container.remove(vnf_k)
                    host.cpu_res += cpu_brc_per
                    host.mem_res += mem_brc_per
                #host_vnf_container_2 = copy.deepcopy(host.vnf_container)
                #print "after restore vnf in host_vnf_container",host_vnf_container_2
            
            return 0
        
    
    def verify_band_left_0(self, path, core_link, NFR, host): #the host does not contain the related NFCR
        
        if NFR[1].pre_vnf_id == None and NFR[1].next_vnf_id == None: #the NFCR contains only one VNF
            host.verify_band(NFR[1].in_band)
        if NFR[1].pre_vnf_id == None and NFR[1].next_vnf_id != None:  #the VNF is the first
            host.verify_band(NFR[1].in_band)
        if NFR[1].pre_vnf_id != None and NFR[1].next_vnf_id == None: #the VNF is the last
            host.verify_band(NFR[1].in_band)
        if NFR[1].pre_vnf_id != None and NFR[1].next_vnf_id != None: #the VNF in the middle
            host.verify_band(NFR[1].in_band + NFR[1].out_band)
    
    def verify_band_left_1(self, path, core_link, NFR, host):  #the host contains the related NFCR
        vnf_in_NFCR = self.calc_req_vnf(host.req_container_dict[NFR[0]])
        if NFR[1].pre_vnf_id in vnf_in_NFCR and NFR[1].next_vnf_id not in vnf_in_NFCR: 
            host.band_res += NFR[1].in_band
            host.verify_band(NFR[1].out_band)
        if NFR[1].pre_vnf_id not in vnf_in_NFCR and NFR[1].next_vnf_id in vnf_in_NFCR:
            host.band_res += NFR[1].out_band
            host.verify_band(NFR[1].in_band)
        if NFR[1].pre_vnf_id not in vnf_in_NFCR and NFR[1].next_vnf_id not in vnf_in_NFCR:
            host.verify_band(NFR[1].in_band + NFR[1].out_band)
        if NFR[1].pre_vnf_id in vnf_in_NFCR and NFR[1].next_vnf_id in vnf_in_NFCR:
            host.band_res += NFR[1].in_band + NFR[1].out_band
    
    def get_aver(self, value_list):
        list_len = len(value_list)
        v_total = 0
        for i in range(list_len):
            v_total += value_list[i]
        aver = v_total/list_len
        return aver
    
    def calc_host_resource(self, host_candinate_list):
        host_resource = {}
        for i in range(len(host_candinate_list)):
            h_value = host_candinate_list[i]
            host_resource[h_value.id] = self.get_aver(h_value.cpu_res)/CPU_TOTAL + self.get_aver(h_value.mem_res)/MEM_TOTAL + self.get_aver(h_value.band_res)/h_band
        #min_h_id = self.search_min_dict(host_resource)
        return host_resource
    
    def search_min_dict(self, dict_exp):
        min_val = INFINITY
        for key, value in dict_exp.iteritems():
            if value < min_val:
                min_val = value
                min_id = key
        return min_id
    
    
    def search_max_dict(self, dict_exp):
        max_val = 0
        for key, value in dict_exp.iteritems():
            if value > max_val:
                max_val = value
                max_id = key
        return max_id
    
    def sort_requests(self, NFR_lib):
        req_aver_source_dict = {}
        for i in range(len(NFR_lib)):
            VNF_i_cpu = 0
            VNF_i_mem = 0
            VNF_i_band = 0
            VNF_i_cpu += self.get_aver(NFR_lib[i][1].cpu)/CPU_TOTAL
            VNF_i_mem += self.get_aver(NFR_lib[i][1].mem)/MEM_TOTAL
            VNF_i_band += self.get_aver(NFR_lib[i][1].in_band + NFR_lib[i][1].out_band)/h_band  #the band resource deserves recosideration
            req_aver_source_dict[i] = VNF_i_cpu + VNF_i_mem + VNF_i_band
        req_id_list = []
        while(len(req_aver_source_dict)):
            max_id = self.search_max_dict(req_aver_source_dict)
            req_id_list.append(max_id)
            req_aver_source_dict.pop(max_id)
        
        return req_id_list
    
    def algorithm_h(self, route, core_link, e_link, request_lib, host_lib, request_id_list):      #algorithm_h
        print "algorithm_h no adjustment"
        
        host_id = 0
        host_list = []
        req_abandoned = []
        host_list.append(host_lib[host_id])
        request_num = len(request_lib)
        sw_id = route._get_switch(host_id,e_link)
        path = route._get_path(ROOT_ID, sw_id)
        #NFR_lib_backup = copy.deepcopy(NFR_lib)
        for k in request_id_list:
            NFR_lib = {}
            NFR_num = 0
            for NFR_i in request_lib[k]:
                NFR_tuple = (k,NFR_i)
                NFR_lib[NFR_num] = NFR_tuple
                NFR_num += 1
            NFR_id_list = self.sort_requests(NFR_lib)
            
            for i in NFR_id_list:
                j = 0
                host_candinate_list = []
                for j in range(len(host_list)):
                    if self.check_end(path, core_link, NFR_lib[i], host_list[j]):
                        host_candinate_list.append(host_list[j])
                if len(host_candinate_list):
                    host_resource = self.calc_host_resource(host_candinate_list)
                    host_chosen_id = self.search_min_dict(host_resource)
                    
                    if NFR_lib[i][0] in host_lib[host_chosen_id].req_container:
                        host_lib[host_chosen_id].req_container_dict[NFR_lib[i][0]].append(NFR_lib[i][1])
                        self.verify_node_resource_left(path, core_link, NFR_lib[i],host_lib[host_chosen_id])
                        self.verify_band_left_1(path,core_link,NFR_lib[i],host_lib[host_chosen_id])
                    else:
                        host_lib[host_chosen_id].req_container.append(NFR_lib[i][0])
                        host_lib[host_chosen_id].req_container_dict[NFR_lib[i][0]] = [NFR_lib[i][1]]
                        self.verify_node_resource_left(path, core_link, NFR_lib[i],host_lib[host_chosen_id])
                        self.verify_band_left_0(path,core_link,NFR_lib[i],host_lib[host_chosen_id])
                
                if len(host_candinate_list) == 0:
                    host_id += 1
                    host_list.append(host_lib[host_id])
                    sw_id = route._get_switch(host_id,e_link)
                    path = route._get_path(ROOT_ID, sw_id)
                    if self.check_end(path, core_link, NFR_lib[i], host_lib[host_id]) == 1:
                        if NFR_lib[i][0] in host_lib[host_id].req_container:
                            host_lib[host_id].req_container_dict[NFR_lib[i][0]].append(NFR_lib[i][1])
                            self.verify_node_resource_left(path, core_link, NFR_lib[i],host_lib[host_id])
                            self.verify_band_left_1(path,core_link,NFR_lib[i],host_lib[host_id])
                        else:
                            host_lib[host_id].req_container.append(NFR_lib[i][0])
                            host_lib[host_id].req_container_dict[NFR_lib[i][0]] = [NFR_lib[i][1]]
                            self.verify_node_resource_left(path, core_link, NFR_lib[i],host_lib[host_id])
                            self.verify_band_left_0(path,core_link,NFR_lib[i],host_lib[host_id])
                    else:
                        req_abandoned.append(i)
        
        print "req_abandoned is:",req_abandoned,'\n',len(req_abandoned)
        host_obsoleted = []
        for k in range(host_id+1):
            if len(host_lib[k].req_container) == 0:
                host_obsoleted.append(k)
        host_id = host_id - len(host_obsoleted)
        print "final host id (Algorithm H no adjustment):", host_id
        return (host_id,host_lib)


class algorithm_h_adjustment(object):
    def search_max_dict(self, dict_exp):
        max_val = 0
        for key, value in dict_exp.iteritems():
            if value > max_val:
                max_val = value
                max_id = key
        return max_id
    
    
    def get_aver(self, value_list):
        list_len = len(value_list)
        v_total = 0
        for i in range(list_len):
            v_total += value_list[i]
        aver = v_total/list_len
        return aver
    
    def search_min_dict(self, dict_exp):
        min_val = INFINITY
        for key, value in dict_exp.iteritems():
            if value < min_val:
                min_val = value
                min_id = key
        return min_id
    
    def calc_host_resource(self, host_id, host_lib):
        host_resource = {}
        for i in range(host_id+1):
            h_value = host_lib[i]
            host_resource[i] = self.get_aver(h_value.cpu_res)/CPU_TOTAL + self.get_aver(h_value.mem_res)/MEM_TOTAL + self.get_aver(h_value.band_res)/h_band
        #min_h_id = self.search_min_dict(host_resource)
        return host_resource
    
    def sort_host(self, host_resource):
        h_sorted_list = []
        while(len(host_resource)):
            h_max_id = self.search_max_dict(host_resource)
            h_sorted_list.append(h_max_id)
            host_resource.pop(h_max_id)
        return h_sorted_list    
     
    
    def verify_resource_left_dst(self, request, host):
        req_len = len(request)
        req_cpu = 0
        req_mem = 0
        req_band = 0
        for i in range(req_len):
            req_cpu += request[i].cpu
            req_mem += request[i].mem
        #req_band = request[0].band                            #request band is determined by the first VNF
        
        host.verify_cpu(req_cpu)
        host.verify_mem(req_mem)
        
        
    def verify_resource_left_src(self, request, host):
        req_len = len(request)
        req_cpu = 0
        req_mem = 0
        req_band = 0
        for i in range(req_len):
            req_cpu += request[i].cpu
            req_mem += request[i].mem
        #req_band = request[0].band                            #request band is determined by the first VNF
        host.cpu_res += req_cpu
        host.mem_res += req_mem
        
  
    def check_value(self, list_ori, list_cmp):
        list_len = len(list_ori)
        break_t = 0
        for i in range(list_len):
            if list_cmp[i] >= list_ori[i]:
                break_t += 1
        return break_t

    def check_host(self, vnf_temp, host_dst):
        req_cpu = vnf_temp.cpu
        req_mem = vnf_temp.mem
        req_band = vnf_temp.in_band + vnf_temp.out_band
        #req_band = vnf_temp.in_band
        
        h_cpu_bt_total = self.check_value(host_dst.cpu_res, req_cpu)                       #break times
        h_mem_bt_total = self.check_value(host_dst.mem_res, req_mem)
        h_band_bt_total = self.check_value(host_dst.band_res, req_band)
        
        if h_cpu_bt_total or h_mem_bt_total or h_band_bt_total:                     #to be changed, a threshold can be defined
            return False
        else:
            host_dst.verify_band(req_band)
            host_dst.verify_cpu(req_cpu)
            host_dst.verify_mem(req_mem)
            return True
    
    def verify_band_and_req_dict(self, transformed_request, host_dst, host_src):
        
        def vnfID_in_req(vnf_list):
            vnfID_list = []
            for vnf_i in vnf_list:
                vnfID_list.append(vnf_i.id)
            return vnfID_list
            
        
        #the host_dst has no the same req
        def vnf_first_0(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id): #the req has not been moved clear
                vnf_id_list = vnfID_in_req(src_host.req_container_dict[tran_id])
                if trans_vnf.next_vnf_id in vnf_id_list:       #next_vnf_id in the src_host
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                else:            #next_vnf_id in the third part, because the host_dst has no the same req_id,so it has no the next_vnf_id definitely.
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
            else:
                #the next_vnf_id is not none, and the next_vnf_id is in the third part,because the host_dst has no the same req_id,so it has no the next_vnf_id definitely.
                dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
        
        def vnf_last_0(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list = vnfID_in_req(src_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list:
                    dst_host.band_res -= trans_vnf.in_band
                    src_host.band_res -= trans_vnf.in_band
                else:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res -= trans_vnf.in_band
            else:
                src_host.band_res += trans_vnf.in_band
                dst_host.band_res -= trans_vnf.in_band
        
        def vnf_middle_0(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list = vnfID_in_req(src_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list and trans_vnf.next_vnf_id in vnf_id_list:
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list and trans_vnf.next_vnf_id not in vnf_id_list:
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id not in vnf_id_list and trans_vnf.next_vnf_id in vnf_id_list:
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                else:
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
            else:
                dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
        
        #the host_dst has the same req
        def vnf_first_1(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list_src = vnfID_in_req(src_host.req_container_dict[tran_id])
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.next_vnf_id in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                else:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
            else:
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                else:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
        
        def vnf_last_1(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list_src = vnfID_in_req(src_host.req_container_dict[tran_id])
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list_src:
                    src_host.band_res -= trans_vnf.in_band
                    dst_host.band_res -= trans_vnf.in_band
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res += trans_vnf.in_band
                else:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res -= trans_vnf.in_band
            else:
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list_dst:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res += trans_vnf.in_band
                else:
                    src_host.band_res += trans_vnf.in_band
                    dst_host.band_res -= trans_vnf.in_band
        
        def vnf_middle_1(dst_host, src_host, tran_id, trans_vnf):
            if src_host.req_container_dict.has_key(tran_id):
                vnf_id_list_src = vnfID_in_req(src_host.req_container_dict[tran_id])
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list_src and trans_vnf.next_vnf_id in vnf_id_list_src:
                    src_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_src and trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_src and trans_vnf.next_vnf_id not in vnf_id_list_dst and trans_vnf.next_vnf_id not in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.out_band - trans_vnf.in_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id not in vnf_id_list_dst and trans_vnf.next_vnf_id not in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.out_band + trans_vnf.in_band)
                    dst_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.next_vnf_id in vnf_id_list_src and trans_vnf.pre_vnf_id not in vnf_id_list_dst and trans_vnf.pre_vnf_id not in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.next_vnf_id in vnf_id_list_dst and trans_vnf.pre_vnf_id not in vnf_id_list_dst and trans_vnf.pre_vnf_id not in vnf_id_list_src:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.out_band - trans_vnf.in_band)
                else:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
            else:
                vnf_id_list_dst = vnfID_in_req(dst_host.req_container_dict[tran_id])
                if trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                elif trans_vnf.pre_vnf_id in vnf_id_list_dst and trans_vnf.next_vnf_id not in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res += (trans_vnf.in_band - trans_vnf.out_band)
                elif trans_vnf.next_vnf_id in vnf_id_list_dst and trans_vnf.pre_vnf_id not in vnf_id_list_dst:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band - trans_vnf.out_band)
                else:
                    src_host.band_res += (trans_vnf.in_band + trans_vnf.out_band)
                    dst_host.band_res -= (trans_vnf.in_band + trans_vnf.out_band)
        
        for trans_id, trans_value in transformed_request.iteritems():
            if host_dst.req_container_dict.has_key(trans_id):
                #self.verify_req_container_src(host_src)
                host_dst.req_container_dict[trans_id].append(trans_value)
                if trans_value.pre_vnf_id == None:
                    vnf_first_1(host_dst, host_src, trans_id, trans_value)
                elif trans_value.next_vnf_id == None:
                    vnf_last_1(host_dst, host_src, trans_id, trans_value)
                else:
                    vnf_middle_1(host_dst, host_src, trans_id, trans_value)
                #host_dst.req_container_dict[trans_id] = self.right_sequence_vnf(host_dst.req_container_dict[trans_id], requests_lib[trans_id])
            else:
                #self.verify_req_container_src(host_src)
                host_dst.req_container_dict[trans_id] = [trans_value]
                host_dst.req_container.append(trans_id)
                if trans_value.pre_vnf_id == None and trans_value.next_vnf_id == None:
                    host_dst.band_res -= trans_value.in_band
                    host_src.band_res += trans_value.in_band
                elif trans_value.pre_vnf_id == None and trans_value.next_vnf_id != None:
                    vnf_first_0(host_dst, host_src, trans_id, trans_value)
                elif trans_value.next_vnf_id == None and trans_value.pre_vnf_id != None:
                    vnf_last_0(host_dst, host_src, trans_id, trans_value)
                else:
                    vnf_middle_0(host_dst, host_src, trans_id, trans_value)
   
    
    def count_species(self, host_lib):           #count the vnf species in the host
        VNF_list = []
        for req_id,req_value in host_lib.req_container_dict.iteritems():
            for VNF in req_value:
                if VNF.id in VNF_list:
                    pass
                else:
                    VNF_list.append(VNF.id)
            
        host_lib.vnf_container = copy.deepcopy(VNF_list)          #each host has a vnf list to hold all the different vnf
        return VNF_list


    def algorithm_h_adjustment(self, route, core_link, e_link, host_id, host_lib):
        host_empty = []
        host_non_empty = []
        for i in range(host_id+1):
            if len(host_lib[i].vnf_container) == 0:                                 
                host_empty.append(i)
            else:
                host_non_empty.append(i)
        print "****start to calc the rest resource in each host****"
        host_resource = self.calc_host_resource(host_id, host_lib)
        print "****start to sort the host list based on aver resource****"
        h_sorted_list = self.sort_host(host_resource)
        hsl_len = len(h_sorted_list)
        
        start_p = len(host_empty)
        end_p = start_p + 1
        #print "****the start_p and end_p are****",start_p, end_p
        shift_p = 0
        abandoned_host_id = []
        host_src_id = h_sorted_list[start_p]
        host_src = host_lib[host_src_id]
        host_lib_backup = copy.deepcopy(host_lib)
        while(start_p < end_p):
            #print "****the start_p and end_p are****",start_p, end_p
            host_dst_id = h_sorted_list[end_p]
            host_dst = host_lib[host_dst_id]
            vnf_len_in_src_ori = self.count_species(host_src)
            
            #print "req_in_host_dst is:",req_in_host_dst
            #sw_id = route._get_switch(host_id,e_link)
            #path = route._get_path(ROOT_ID, sw_id)
            
            req_removed_list = []
            host_dst_for_check = copy.deepcopy(host_dst)
            for req_src_id,req_src_value in host_src.req_container_dict.iteritems():
                for vnf_index in req_src_value:
                    if vnf_index.id in host_dst.vnf_container and self.check_host(vnf_index,host_dst_for_check):
                        req_removed_list.append((req_src_id,vnf_index))
                        req_list = [vnf_index]
                        self.verify_resource_left_src(req_list, host_src)
                        self.verify_resource_left_dst(req_list, host_dst)
            
            for req_id_removed in req_removed_list:
                trans_req = {}
                trans_req[req_id_removed[0]] = req_id_removed[1]
                    
                for vnf_index in host_src.req_container_dict[req_id_removed[0]]:
                    if vnf_index.id == req_id_removed[1].id:
                        vnf_removed = vnf_index
                    
                host_src.req_container_dict[req_id_removed[0]].remove(vnf_removed)
                if len(host_src.req_container_dict[req_id_removed[0]]) == 0:
                    host_src.req_container_dict.pop(req_id_removed[0])
                    host_src.req_container.remove(req_id_removed[0])
                self.verify_band_and_req_dict(trans_req, host_dst, host_src)
            vnf_len_in_src_new = self.count_species(host_src)
            host_src.cpu_res += (len(vnf_len_in_src_ori) - len(vnf_len_in_src_new))*cpu_brc_per
            host_src.mem_res += (len(vnf_len_in_src_ori) - len(vnf_len_in_src_new))*mem_brc_per
            for vnf_id_ori in vnf_len_in_src_ori:
                if vnf_id_ori in vnf_len_in_src_new:
                    continue
                else:
                    host_src.vnf_pop_container.append(vnf_id_ori)
            
            if len(host_src.req_container_dict.keys()) == 0:
                print "host_src", host_src_id, "has been moved clear"
                host_empty.append(host_src_id)
                print "****the start_p and end_p are(len(host_src.req_container) == 0)****",start_p, end_p
                break
            else:
                end_p -= 1
        if start_p >= end_p:
            host_lib = copy.deepcopy(host_lib_backup)
            print "saved host after adjustment is:",len(host_empty)
            return (host_empty,host_lib)
        else:
            algorithm_h_adjustment(route, core_link, e_link, host_id, host_lib)

    


def Algorithm_H_r(data_back):
    fattree = fat_tree()
    fattree.K = 16
    #fattree.__init__()
    fattree.instant_topo()
    sws = fattree.switch
    adj = fattree.core_link
    
    
    core_link_AH = copy.deepcopy(fattree.core_link)
    e_link_AH = copy.deepcopy(fattree.edge_link)

    
    
    request_lib_AH = copy.deepcopy(data_back)
    host_lib_AH = copy.deepcopy(fattree.host_dict)
    
    route = VNF_route()
    route.switches = sws
    route.adjacency = adj
    #print "*******test route*******"
    #route._get_path(7,18)
    #print "test route over"
    req_2_host = VNF_algorithm_h()
    req_2_host_adjust = algorithm_h_adjustment()
    #print "******pre process of request_lib in order to reset all the value in the list as peak value*****"
    #pre_process_req_lib(request_lib)
    #print "*******pre process is over******"
    print "*****start to map the requests to hosts*****"
    req_sequence_list = []
    for i in range(len(request_lib_AH)):
        req_sequence_list.append(i)
    
    random.shuffle(req_sequence_list)
    #host_id = 128
    host_info_AH = req_2_host.algorithm_h(route, core_link_AH, e_link_AH, request_lib_AH, host_lib_AH, req_sequence_list)
    
    host_info_AH_AD = req_2_host_adjust.algorithm_h_adjustment(route, core_link_AH, e_link_AH, host_info_AH[0], host_info_AH[1])
    
    print "*********start save info*********"
    host_id_AH = host_info_AH[0] - len(host_info_AH_AD[0])
    host_lib_AH = host_info_AH_AD[1]
    
    print "the used host num(AH) before adjustment is:",host_info_AH[0]+1
    print "the used host num(AH) after adjustment is:",host_id_AH+1
    
    #with open('G:/ldf/simulation/data_store/host_lib_AH_r.pickle', 'w') as f:
        #pickle.dump(host_lib_AH,f)
    saved_info_dict_AH = {}
    saved_info_dict_AH[0] = host_id_AH+1
    with open('G:/ldf/simulation/data_store/saved_info_dict_AH_r.json', 'w') as f:
        json.dump(saved_info_dict_AH,f)
    
    
    print "***********Algorithm_H_r is over**********"
 
 
 
def Algorithm_H_nr(data_back):
    fattree = fat_tree()
    fattree.K = 16
    #fattree.__init__()
    fattree.instant_topo()
    sws = fattree.switch
    adj = fattree.core_link
    
    
    core_link_AH = copy.deepcopy(fattree.core_link)
    e_link_AH = copy.deepcopy(fattree.edge_link)

    
    
    request_lib_AH = copy.deepcopy(data_back)
    host_lib_AH = copy.deepcopy(fattree.host_dict)
    
    route = VNF_route()
    route.switches = sws
    route.adjacency = adj
    #print "*******test route*******"
    #route._get_path(7,18)
    #print "test route over"
    req_2_host = VNF_algorithm_h()
    req_2_host_adjust = algorithm_h_adjustment()
    #print "******pre process of request_lib in order to reset all the value in the list as peak value*****"
    #pre_process_req_lib(request_lib)
    #print "*******pre process is over******"
    print "*****start to map the requests to hosts*****"
    req_sequence_list = []
    for i in range(len(request_lib_AH)):
        req_sequence_list.append(i)
    
    #random.shuffle(req_sequence_list)
    #host_id = 128
    host_info_AH = req_2_host.algorithm_h(route, core_link_AH, e_link_AH, request_lib_AH, host_lib_AH, req_sequence_list)
    
    host_info_AH_AD = req_2_host_adjust.algorithm_h_adjustment(route, core_link_AH, e_link_AH, host_info_AH[0], host_info_AH[1])
    
    print "*********start save info*********"
    host_id_AH = host_info_AH[0] - len(host_info_AH_AD[0])
    host_lib_AH = host_info_AH_AD[1]
    
    print "the used host num(AH) before adjustment is:",host_info_AH[0]+1
    print "the used host num(AH) after adjustment is:",host_id_AH+1
    
    #with open('G:/ldf/simulation/data_store/host_lib_AH_nr.pickle', 'w') as f:
        #pickle.dump(host_lib_AH,f)
    saved_info_dict_AH = {}
    saved_info_dict_AH[0] = host_id_AH+1
    with open('G:/ldf/simulation/data_store/saved_info_dict_AH_nr.json', 'w') as f:
        json.dump(saved_info_dict_AH,f)
    
    
    print "***********Algorithm_H_nr is over**********"
 
 