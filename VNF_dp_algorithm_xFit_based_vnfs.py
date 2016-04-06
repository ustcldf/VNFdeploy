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

#x = np.arange(0,24,0.1)
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
            print
      
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
    
  
class VNF_xFit(object):
    
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
    
    def NF_map(self, route, core_link, e_link, NFR_lib, host_lib, NFR_id_list):     #map the request to host using Next Fit
        print "Next Fit"
        host_id = 0
        req_abandoned = []
        request_num = len(NFR_lib)
        #sw_id = route._get_switch(host_id,e_link)
        #path = route._get_path(ROOT_ID, sw_id)
        for i in NFR_id_list:
            sw_id = route._get_switch(host_id,e_link)
            path = route._get_path(ROOT_ID, sw_id)
            if self.check_end(path, core_link, NFR_lib[i], host_lib[host_id]):
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
                host_id += 1
                sw_id = route._get_switch(host_id,e_link)
                path = route._get_path(ROOT_ID, sw_id)
                if self.check_end(path, core_link, NFR_lib[i], host_lib[host_id]):
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
        print "final host id (NF):", host_id
        print "abandoned requeats are:",req_abandoned
        return (host_id,host_lib)
    
    def FF_map(self, route, core_link, e_link, NFR_lib, host_lib, NFR_id_list):      #map the request to host using First Fit
        print "First Fit"
        req_abandoned = []
        host_id = 0
        host_list = []
        host_list.append(host_lib[host_id])
        request_num = len(NFR_lib)
        sw_id = route._get_switch(host_id,e_link)
        path = route._get_path(ROOT_ID, sw_id)
        #NFR_lib_backup = copy.deepcopy(NFR_lib)
        for i in NFR_id_list:
            j = 0
            inditor = 0
            for j in range(len(host_list)):
                #print "request seq_num:", i, '\t', "host_id now", host_list[j].id
                if self.check_end(path, core_link, NFR_lib[i], host_list[j]) == 1:
                    if NFR_lib[i][0] in host_lib[host_list[j].id].req_container:
                        host_lib[host_list[j].id].req_container_dict[NFR_lib[i][0]].append(NFR_lib[i][1])
                        self.verify_node_resource_left(path, core_link, NFR_lib[i],host_list[j])
                        self.verify_band_left_1(path,core_link,NFR_lib[i],host_list[j])
                    else:
                        host_lib[host_list[j].id].req_container.append(NFR_lib[i][0])
                        host_lib[host_list[j].id].req_container_dict[NFR_lib[i][0]] = [NFR_lib[i][1]]
                        self.verify_node_resource_left(path, core_link, NFR_lib[i],host_list[j])
                        self.verify_band_left_0(path,core_link,NFR_lib[i],host_list[j])
                    inditor = 1                  #an inditor, indite if there is a need to start a new host
                    break                       #find the first host to hold the request, and stop the process
            if inditor == 0:
                host_id += 1
                host_list.append(host_lib[host_id])
                #print "request seq_num:", i, '\t', "host_id now", host_id 
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
        print "final host id (FF):", host_id
        return (host_id,host_lib)
                
            
        
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
    
    def BF_map(self, route, core_link, e_link, NFR_lib, host_lib, NFR_id_list):      #map the request to host using Best Fit
        print "Best Fit"
        
        host_id = 0
        host_list = []
        req_abandoned = []
        host_list.append(host_lib[host_id])
        request_num = len(NFR_lib)
        sw_id = route._get_switch(host_id,e_link)
        path = route._get_path(ROOT_ID, sw_id)
        #NFR_lib_backup = copy.deepcopy(NFR_lib)
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
        print "final host id (BF):", host_id
        return (host_id,host_lib)
    
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
    
    def FFD_map(self, route, core_link, e_link, request_lib, host_lib):
        print "FFD"
        print "sort the requests based on the aver resource"
        req_id_list = self.sort_requests(request_lib)
        host_info = self.FF_map(route, core_link, e_link, request_lib, host_lib, req_id_list)
        return host_info
    
    def BFD_map(self, route, core_link, e_link, request_lib, host_lib):
        print "BFD"
        print "sort the requests based on the aver resource"
        req_id_list = self.sort_requests(request_lib)
        host_info = self.BF_map(route, core_link, e_link, request_lib, host_lib, req_id_list)
        return host_info
    
    
 


def NFCR_to_NFRs(request_lib):
    NFR_num = 0
    NFR_dict = {}
    for NFCR_id,NFCR_value in request_lib.iteritems():
        for NFR_i in NFCR_value:
            NFR_tuple = (NFCR_id,NFR_i)
            NFR_dict[NFR_num] = NFR_tuple
            NFR_num += 1
    return NFR_dict


    

    
def search_list_max(list_temp):
    max_value = 0
    for i in range(len(list_temp)):
        if list_temp[i] > max_value:
            max_value = list_temp[i]
    return max_value

def search_list_min(list_temp):
    min_value = INFINITY
    for i in range(len(list_temp)):
        if list_temp[i] < min_value:
            min_value = list_temp[i]
    return min_value

def reset_list(list_temp, new_value):
    for i in range(len(list_temp)):
        list_temp[i] = new_value


def pre_process_req_lib(request_lib):
    for req_id,req_value in request_lib.iteritems():
        max_cpu = 0
        max_band = 0
        max_mem = 0
        for i in range(len(req_value)):
            max_cpu = search_list_max(req_value[i].cpu)
            max_mem = search_list_max(req_value[i].mem)
            max_band = search_list_max(req_value[i].band)
            reset_list(req_value[i].cpu, max_cpu)
            reset_list(req_value[i].mem, max_mem)
            reset_list(req_value[i].band, max_band)
    
def check_resource_lefted_list(list_temp):
    verified_times = 0
    for i in range(len(list_temp)):
        if list_temp[i] <= 0:
            verified_times += 1
    return verified_times

def get_aver_list(list_temp):
    aver_value = 0
    total_value = 0
    for i in range(len(list_temp)):
        total_value += list_temp[i]
    aver_value = total_value/len(list_temp)
    return aver_value
    



def after_pod_adjust_check(host_lib, host_id_ori):
    #host_empty_list = []
    print "let's start to check the data in the host_lib"
    print "ID - host_id, A - len of req_container, B - keys of req_container_dict, C - len of vnf_container"
    print "D1 - max of cpu lefted + cpu_in_req, D2 - min of cpu lefted + cpu_in_req, D3 - aver of cpu lefted + cpu_in_req"
    print "E1 - max of mem lefted + mem_in_req, E2 - min of mem lefted + mem_in_req, E3 - aver of mem lefted + mem_in_req"
    print "F1 - max of band lefted + band_in_req, F2 - min of band lefted + band_in_req, F3 - aver of band lefted + band_in_req"
    print "ID**A**B**C**D1**D2**D3**E1**E2**E3**F1**F2**F3"
    for i in range(host_id_ori+1):
        #cpu_verified_t = check_resource_lefted_list(host_lib[i].cpu_res)
        #mem_verified_t = check_resource_lefted_list(host_lib[i].mem_res)
        #band_verified_t = check_resource_lefted_list(host_lib[i].band_res)
        A = len(host_lib[i].req_container)
        #B = len(host_lib[i].req_container_dict)
        B = host_lib[i].req_container_dict.keys()
        #C = len(host_lib[i].vnf_container)
        C = host_lib[i].vnf_container
        req_cpu = 0
        req_mem = 0
        req_band = 0
        if len(host_lib[i].req_container_dict.keys()):
            for req_id, req_value in host_lib[i].req_container_dict.iteritems():
                if len(req_value):
                    for vnf in req_value:
                        req_cpu += vnf.cpu
                        req_mem += vnf.mem
                    req_band += req_value[0].in_band
          
        
        
        D3 = cpu_aver_left = get_aver_list(host_lib[i].cpu_res + req_cpu + len(host_lib[i].vnf_container) * cpu_brc_per)
        E3 = mem_aver_left = get_aver_list(host_lib[i].mem_res + req_mem + len(host_lib[i].vnf_container) * mem_brc_per)
        F3 = band_aver_left = get_aver_list(host_lib[i].band_res + req_band)
        D1 = cpu_max_left = search_list_max(host_lib[i].cpu_res + req_cpu + len(host_lib[i].vnf_container) * cpu_brc_per)
        E1 = mem_max_left = search_list_max(host_lib[i].mem_res + req_mem + len(host_lib[i].vnf_container) * mem_brc_per)
        F1 = band_max_left = search_list_max(host_lib[i].band_res + req_band)
        D2 = cpu_min_left = search_list_min(host_lib[i].cpu_res + req_cpu + len(host_lib[i].vnf_container) * cpu_brc_per)
        E2 = mem_min_left = search_list_min(host_lib[i].mem_res + req_mem + len(host_lib[i].vnf_container) * mem_brc_per)
        F2 = band_min_left = search_list_min(host_lib[i].band_res + req_band)
        
        print i,'**',A,'**',B,'**',C,'**',D1,'**',D2,'**',D3,'**',E1,'**',E2,'**',E3,'**',F1,'**',F2,'**',F3


def host_resource_left_check(host_lib, host_id_ori, host_saved_list):
    host_empty_list = []
    print "let's start to check the data in the host_lib"
    print "ID - host_id, A - len of req_container, B - keys of req_container_dict, C - len of vnf_container"
    print "D1 - max of cpu lefted, D2 - min of cpu lefted, D3 - aver of cpu lefted"
    print "E1 - max of mem lefted, E2 - min of mem lefted, E3 - aver of mem lefted"
    print "F1 - max of band lefted, F2 - min of band lefted, F3 - aver of band lefted"
    print "ID**A**B**C**D1**D2**D3**E1**E2**E3**F1**F2**F3"
    for i in range(host_id_ori+1):
        #cpu_verified_t = check_resource_lefted_list(host_lib[i].cpu_res)
        #mem_verified_t = check_resource_lefted_list(host_lib[i].mem_res)
        #band_verified_t = check_resource_lefted_list(host_lib[i].band_res)
        A = len(host_lib[i].req_container)
        B = len(host_lib[i].req_container_dict)
        #B = host_lib[i].req_container_dict.keys()
        C = len(host_lib[i].vnf_container)
        #C = host_lib[i].vnf_container
        D3 = cpu_aver_left = get_aver_list(host_lib[i].cpu_res)
        E3 = mem_aver_left = get_aver_list(host_lib[i].mem_res)
        F3 = band_aver_left = get_aver_list(host_lib[i].band_res)
        D1 = cpu_max_left = search_list_max(host_lib[i].cpu_res)
        E1 = mem_max_left = search_list_max(host_lib[i].mem_res)
        F1 = band_max_left = search_list_max(host_lib[i].band_res)
        D2 = cpu_min_left = search_list_min(host_lib[i].cpu_res)
        E2 = mem_min_left = search_list_min(host_lib[i].mem_res)
        F2 = band_min_left = search_list_min(host_lib[i].band_res)
        
        print i,'**',A,'**',B,'**',C,'**',D1,'**',D2,'**',D3,'**',E1,'**',E2,'**',E3,'**',F1,'**',F2,'**',F3
        if len(host_lib[i].req_container) == 0:
            host_empty_list.append(i)
        #print "the host has been checked is:",i
        #print "the verified times of cpu", cpu_verified_t
        #print "the verified times of mem", mem_verified_t
        #print "the verified times of band", band_verified_t
        #print "the aver value of lefted cpu", cpu_aver_left
        #print "the aver value of lefted mem", mem_aver_left
        #print "the aver value of lefted band", band_aver_left
    """
    print "Now, we compare the empty list"
    for host_id in host_empty_list:
        if host_id in host_saved_list:
            print host_id, "is in host_saved_list"
        else:
            print "the host_id", host_id,"is not in host_saved_list"
            print "the host_empty_list is not the same with host_saved_list"
    """
    

def VNF_simulation_xFit_randomized_w(data_back):
    fattree = fat_tree()
    fattree.K = 16
    #fattree.__init__()
    fattree.instant_topo()
    sws = fattree.switch
    adj = fattree.core_link
    
    core_link_ff = copy.deepcopy(fattree.core_link)
    e_link_ff = copy.deepcopy(fattree.edge_link)
    core_link_ffd = copy.deepcopy(fattree.core_link)
    e_link_ffd = copy.deepcopy(fattree.edge_link)
    core_link_bf = copy.deepcopy(fattree.core_link)
    e_link_bf = copy.deepcopy(fattree.edge_link)
    core_link_bfd = copy.deepcopy(fattree.core_link)
    e_link_bfd = copy.deepcopy(fattree.edge_link)
    core_link_nf = copy.deepcopy(fattree.core_link)
    e_link_nf = copy.deepcopy(fattree.edge_link)
    
    request_lib_ff = copy.deepcopy(data_back)
    host_lib_ff = copy.deepcopy(fattree.host_dict)
    request_lib_ffd = copy.deepcopy(data_back)
    host_lib_ffd = copy.deepcopy(fattree.host_dict)
    request_lib_bf = copy.deepcopy(data_back)
    host_lib_bf = copy.deepcopy(fattree.host_dict)
    request_lib_bfd = copy.deepcopy(data_back)
    host_lib_bfd = copy.deepcopy(fattree.host_dict)
    request_lib_nf = copy.deepcopy(data_back)
    host_lib_nf = copy.deepcopy(fattree.host_dict)
    
    NFR_lib_ff = NFCR_to_NFRs(request_lib_ff)
    NFR_lib_ffd = NFCR_to_NFRs(request_lib_ffd)
    NFR_lib_bf = NFCR_to_NFRs(request_lib_bf)
    NFR_lib_bfd = NFCR_to_NFRs(request_lib_bfd)
    NFR_lib_nf = NFCR_to_NFRs(request_lib_nf)
    
    route = VNF_route()
    route.switches = sws
    route.adjacency = adj
    #print "*******test route*******"
    #route._get_path(7,18)
    #print "test route over"
    req_2_host = VNF_xFit()
    #print "******pre process of request_lib in order to reset all the value in the list as peak value*****"
    #pre_process_req_lib(request_lib)
    #print "*******pre process is over******"
    print "*****start to map the requests to hosts*****"
    print "Original xFit"
    req_sequence_list = []
    for i in range(len(NFR_lib_nf)):
        req_sequence_list.append(i)
    
    random.shuffle(req_sequence_list)
    #host_id = 128
    host_info_nf = req_2_host.NF_map(route, core_link_nf, e_link_nf, NFR_lib_nf, host_lib_nf, req_sequence_list)
    host_info_ff = req_2_host.FF_map(route, core_link_ff, e_link_ff, NFR_lib_ff, host_lib_ff, req_sequence_list)
    host_info_bf = req_2_host.BF_map(route, core_link_bf, e_link_bf, NFR_lib_bf, host_lib_bf, req_sequence_list)
    print "Descending xFit"
    host_info_ffd = req_2_host.FFD_map(route, core_link_ffd, e_link_ffd, NFR_lib_ffd, host_lib_ffd)
    host_info_bfd = req_2_host.BFD_map(route, core_link_bfd, e_link_bfd, NFR_lib_bfd, host_lib_bfd)
    
    print "*********start save info*********"
    
    host_id_nf = host_info_nf[0]
    host_lib_nf = host_info_nf[1]
    
    #with open('G:/ldf/simulation/data_store/host_lib_NF_rw.pickle', 'w') as f:
        #pickle.dump(host_lib_nf,f)
    saved_info_dict_nf = {}
    saved_info_dict_nf[0] = host_id_nf+1
    print "the used host num(nf) is:",host_id_nf+1
    with open('G:/ldf/simulation/data_store/saved_info_dict_NF_rw.json', 'w') as f:
        json.dump(saved_info_dict_nf,f)
    
    #print "check host_lib_nf"
    #saved_host_nf = []
    #host_resource_left_check(host_lib_nf,host_id_nf,saved_host_nf)
    
    host_id_ff = host_info_ff[0]
    host_lib_ff = host_info_ff[1]
    
    #with open('G:/ldf/simulation/data_store/host_lib_FF_rw.pickle', 'w') as f:
        #pickle.dump(host_lib_ff,f)
    saved_info_dict_ff = {}
    saved_info_dict_ff[0] = host_id_ff+1
    print "the used host num(ff) is:",host_id_ff+1
    with open('G:/ldf/simulation/data_store/saved_info_dict_FF_rw.json', 'w') as f:
        json.dump(saved_info_dict_ff,f)
    
    #print "check host_lib_ff"
    #saved_host_ff = []
    #host_resource_left_check(host_lib_ff,host_id_ff,saved_host_ff)
    
    host_id_ffd = host_info_ffd[0]
    host_lib_ffd = host_info_ffd[1]
    print "the used host num(ffd) is:",host_id_ffd+1
    #with open('G:/ldf/simulation/data_store/host_lib_FFD_rw.pickle', 'w') as f:
        #pickle.dump(host_lib_ffd,f)
    saved_info_dict_ffd = {}
    saved_info_dict_ffd[0] = host_id_ffd+1
    with open('G:/ldf/simulation/data_store/saved_info_dict_FFD_rw.json', 'w') as f:
        json.dump(saved_info_dict_ffd,f)
    
    #print "check host_lib_ffd"
    #saved_host_ffd = []
    #host_resource_left_check(host_lib_ffd,host_id_ffd,saved_host_ffd)
    
    host_id_bf = host_info_bf[0]
    host_lib_bf = host_info_bf[1]
    print "the used host num(bf) is:",host_id_bf+1
    
    #with open('G:/ldf/simulation/data_store/host_lib_BF_rw.pickle', 'w') as f:
        #pickle.dump(host_lib_bf,f)
    saved_info_dict_bf = {}
    saved_info_dict_bf[0] = host_id_bf+1
    with open('G:/ldf/simulation/data_store/saved_info_dict_BF_rw.json', 'w') as f:
        json.dump(saved_info_dict_bf,f)
    
    #print "check host_lib_bf"
    #saved_host_bfd = []
    #host_resource_left_check(host_lib_bf,host_id_bf,saved_host_bfd)
    
    host_id_bfd = host_info_bfd[0]
    host_lib_bfd = host_info_bfd[1]
    print "the used host num(bfd) is:",host_id_bfd+1
    
    #with open('G:/ldf/simulation/data_store/host_lib_BFD_rw.pickle', 'w') as f:
        #pickle.dump(host_lib_bfd,f)
    saved_info_dict_bfd = {}
    saved_info_dict_bfd[0] = host_id_bfd+1
    with open('G:/ldf/simulation/data_store/saved_info_dict_BFD_rw.json', 'w') as f:
        json.dump(saved_info_dict_bfd,f)
    
    #print "check host_lib_bfd"
    #saved_host_bfd = []
    #host_resource_left_check(host_lib_bfd,host_id_bfd,saved_host_bfd)
    
    print "***********xFit_random is over**********"
    
def VNF_simulation_xFit_non_randomized_w(data_back):
    fattree = fat_tree()
    fattree.K = 16
    #fattree.__init__()
    fattree.instant_topo()
    sws = fattree.switch
    adj = fattree.core_link
    
    core_link_ff = copy.deepcopy(fattree.core_link)
    e_link_ff = copy.deepcopy(fattree.edge_link)
    core_link_ffd = copy.deepcopy(fattree.core_link)
    e_link_ffd = copy.deepcopy(fattree.edge_link)
    core_link_bf = copy.deepcopy(fattree.core_link)
    e_link_bf = copy.deepcopy(fattree.edge_link)
    core_link_bfd = copy.deepcopy(fattree.core_link)
    e_link_bfd = copy.deepcopy(fattree.edge_link)
    core_link_nf = copy.deepcopy(fattree.core_link)
    e_link_nf = copy.deepcopy(fattree.edge_link)
    
    request_lib_ff = copy.deepcopy(data_back)
    host_lib_ff = copy.deepcopy(fattree.host_dict)
    request_lib_ffd = copy.deepcopy(data_back)
    host_lib_ffd = copy.deepcopy(fattree.host_dict)
    request_lib_bf = copy.deepcopy(data_back)
    host_lib_bf = copy.deepcopy(fattree.host_dict)
    request_lib_bfd = copy.deepcopy(data_back)
    host_lib_bfd = copy.deepcopy(fattree.host_dict)
    request_lib_nf = copy.deepcopy(data_back)
    host_lib_nf = copy.deepcopy(fattree.host_dict)
    
    NFR_lib_ff = NFCR_to_NFRs(request_lib_ff)
    NFR_lib_ffd = NFCR_to_NFRs(request_lib_ffd)
    NFR_lib_bf = NFCR_to_NFRs(request_lib_bf)
    NFR_lib_bfd = NFCR_to_NFRs(request_lib_bfd)
    NFR_lib_nf = NFCR_to_NFRs(request_lib_nf)
    
    
    route = VNF_route()
    route.switches = sws
    route.adjacency = adj
    #print "*******test route*******"
    #route._get_path(7,18)
    #print "test route over"
    req_2_host = VNF_xFit()
    #print "******pre process of request_lib in order to reset all the value in the list as peak value*****"
    #pre_process_req_lib(request_lib)
    #print "*******pre process is over******"
    print "*****start to map the requests to hosts*****"
    print "Original xFit"
    """
    req_sequence_list = [0,1,2,3,4,5,6]  #for test
    for req_id in req_sequence_list:
        vnf_in_req_list = []
        for vnf in request_lib[req_id]:
            vnf_in_req_list.append(vnf.id)
        print "request",req_id,'\t',"vnf num:", vnf_in_req_list
    
    """
    req_sequence_list = []
    for i in range(len(NFR_lib_nf)):
        req_sequence_list.append(i)
    
    #random.shuffle(req_sequence_list)
    #host_id = 128
    host_info_nf = req_2_host.NF_map(route, core_link_nf, e_link_nf, NFR_lib_nf, host_lib_nf, req_sequence_list)
    host_info_ff = req_2_host.FF_map(route, core_link_ff, e_link_ff, NFR_lib_ff, host_lib_ff, req_sequence_list)
    host_info_bf = req_2_host.BF_map(route, core_link_bf, e_link_bf, NFR_lib_bf, host_lib_bf, req_sequence_list)
    print "Descending xFit"
    host_info_ffd = req_2_host.FFD_map(route, core_link_ffd, e_link_ffd, NFR_lib_ffd, host_lib_ffd)
    host_info_bfd = req_2_host.BFD_map(route, core_link_bfd, e_link_bfd, NFR_lib_bfd, host_lib_bfd)

    
    
    print "*********start save info*********"
    host_id_nf = host_info_nf[0]
    host_lib_nf = host_info_nf[1]
    
    with open('G:/ldf/simulation/data_store/host_lib_NF_nrw.pickle', 'w') as f:
        pickle.dump(host_lib_nf,f)
    saved_info_dict_nf = {}
    saved_info_dict_nf[0] = host_id_nf+1
    print "the used host num(nf) is:",host_id_nf+1
    with open('G:/ldf/simulation/data_store/saved_info_dict_NF_nrw.json', 'w') as f:
        json.dump(saved_info_dict_nf,f)
    
    host_id_ff = host_info_ff[0]
    host_lib_ff = host_info_ff[1]
    
    with open('G:/ldf/simulation/data_store/host_lib_FF_nrw.pickle', 'w') as f:
        pickle.dump(host_lib_ff,f)
    saved_info_dict_ff = {}
    saved_info_dict_ff[0] = host_id_ff+1
    print "the used host num(ff) is:",host_id_ff+1
    with open('G:/ldf/simulation/data_store/saved_info_dict_FF_nrw.json', 'w') as f:
        json.dump(saved_info_dict_ff,f)
    
    host_id_ffd = host_info_ffd[0]
    host_lib_ffd = host_info_ffd[1]
    print "the used host num(ffd) is:",host_id_ffd+1
    with open('G:/ldf/simulation/data_store/host_lib_FFD_nrw.pickle', 'w') as f:
        pickle.dump(host_lib_ffd,f)
    saved_info_dict_ffd = {}
    saved_info_dict_ffd[0] = host_id_ffd+1
    with open('G:/ldf/simulation/data_store/saved_info_dict_FFD_nrw.json', 'w') as f:
        json.dump(saved_info_dict_ffd,f)
    
    host_id_bf = host_info_bf[0]
    host_lib_bf = host_info_bf[1]
    print "the used host num(bf) is:",host_id_bf+1
    
    with open('G:/ldf/simulation/data_store/host_lib_BF_nrw.pickle', 'w') as f:
        pickle.dump(host_lib_bf,f)
    saved_info_dict_bf = {}
    saved_info_dict_bf[0] = host_id_bf+1
    with open('G:/ldf/simulation/data_store/saved_info_dict_BF_nrw.json', 'w') as f:
        json.dump(saved_info_dict_bf,f)
    
    host_id_bfd = host_info_bfd[0]
    host_lib_bfd = host_info_bfd[1]
    print "the used host num(bfd) is:",host_id_bfd+1
    
    with open('G:/ldf/simulation/data_store/host_lib_BFD_nrw.pickle', 'w') as f:
        pickle.dump(host_lib_bfd,f)
    saved_info_dict_bfd = {}
    saved_info_dict_bfd[0] = host_id_bfd+1
    with open('G:/ldf/simulation/data_store/saved_info_dict_BFD_nrw.json', 'w') as f:
        json.dump(saved_info_dict_bfd,f)
    #saved_host_bfd = []
    #host_resource_left_check(host_lib_bfd,host_id_bfd,saved_host_bfd)
    print "***********xFit_non_random is over**********"

