#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 14:29:21 2024

@author: li
"""

import copy
import networkx as nx
import random




def purge_smaller_components(network, minimum_number_of_nodes_permitted): #minimum must be int and not lower than 2
    sorted_components = sorted(nx.connected_components(network), key=len, reverse=True)
    for i in range(0, len(sorted_components)):
        if len(sorted_components[i])<minimum_number_of_nodes_permitted:
            purged_network=network.subgraph(set().union(*sorted_components[0:i]))
            break
    return purged_network


def random_disruption_absolute_efficiency(unfrozen_network, number_of_iterations, filename):
    with open(filename, 'a') as f:
        f.write(str(unfrozen_network.name)+',')
    for i in range(0, number_of_iterations):
        if i==number_of_iterations-1:
            print(nx.global_efficiency(unfrozen_network))
            with open(filename, 'a') as f:
                f.write((str(nx.global_efficiency(unfrozen_network)))+'\n')
        else:
            print(nx.global_efficiency(unfrozen_network))
            with open(filename, 'a') as f:
                f.write((str(nx.global_efficiency(unfrozen_network)))+',')
            unfrozen_network.remove_node(random.choice(list(unfrozen_network)))
                
                
def random_disruption_relative_efficiency(unfrozen_network, number_of_iterations, filename):
    with open(filename, 'a') as f:
        f.write(str(unfrozen_network.name)+',')
    original_efficiency=nx.global_efficiency(unfrozen_network)
    for i in range(0, number_of_iterations):
        if i==number_of_iterations-1:
            print(nx.global_efficiency(unfrozen_network)/original_efficiency)
            with open(filename, 'a') as f:
                f.write((str(nx.global_efficiency(unfrozen_network)/original_efficiency))+'\n')
        else:
            print(nx.global_efficiency(unfrozen_network)/original_efficiency)
            with open(filename, 'a') as f:
                f.write((str(nx.global_efficiency(unfrozen_network)/original_efficiency))+',')
            unfrozen_network.remove_node(random.choice(list(unfrozen_network)))
                
                
def random_disruption_absolute_gc_dimension(unfrozen_network, number_of_iterations, filename):
    with open(filename, 'a') as f:
        f.write(str(unfrozen_network.name)+',')
    for i in range(0, number_of_iterations):
        if i==number_of_iterations-1:
            print(len(sorted(nx.connected_components(unfrozen_network), key=len, reverse=True)[0]))
            with open(filename, 'a') as f:
                f.write((str(len(sorted(nx.connected_components(unfrozen_network), key=len, reverse=True)[0]))+'\n'))
        else:
            print(len(sorted(nx.connected_components(unfrozen_network), key=len, reverse=True)[0]))
            with open(filename, 'a') as f:
                f.write((str(len(sorted(nx.connected_components(unfrozen_network), key=len, reverse=True)[0]))+','))
            unfrozen_network.remove_node(random.choice(list(unfrozen_network)))
                
                
def random_disruption_relative_gc_dimension(unfrozen_network, number_of_iterations, filename):
    with open(filename, 'a') as f:
        f.write(str(unfrozen_network.name)+',')
    original_gc_dimension=float(len(sorted(nx.connected_components(unfrozen_network), key=len, reverse=True)[0]))
    for i in range(0, number_of_iterations):
        if i==number_of_iterations-1:
            print(float(len(sorted(nx.connected_components(unfrozen_network), key=len, reverse=True)[0]))/original_gc_dimension)
            with open(filename, 'a') as f:
                f.write((str(float(len(sorted(nx.connected_components(unfrozen_network), key=len, reverse=True)[0]))/original_gc_dimension))+'\n')
        else:
            print(float(len(sorted(nx.connected_components(unfrozen_network), key=len, reverse=True)[0]))/original_gc_dimension)
            with open(filename, 'a') as f:
                f.write((str(float(len(sorted(nx.connected_components(unfrozen_network), key=len, reverse=True)[0]))/original_gc_dimension))+',')
            unfrozen_network.remove_node(random.choice(list(unfrozen_network)))                
                

MC=nx.read_graphml("mcnet_0.65.graphml").to_undirected()
MC.name='MC'
MS=nx.read_graphml("msnet 0.65.graphml").to_undirected()
MS.name='MS'
SC=nx.read_graphml("scnet 0.65.graphml").to_undirected()
SC.name='SC'
SS=nx.read_graphml("ssnet 0.65 0.04.graphml").to_undirected()
SS.name='SS'
HC=nx.read_graphml("hcnet 0.65 0.04.graphml").to_undirected()
HC.name='HC'

normal_networks=[HC, MS, MC, SS, SC]
networks_purgati=[]


for network in normal_networks:
    nx.relabel_nodes(network, nx.get_node_attributes(network, 'shared name'), copy=False)
    network_purgato=purge_smaller_components(network, 5)
    network_purgato_scongelato=nx.Graph(network_purgato)
    networks_purgati.append(network_purgato_scongelato)

    
for network in networks_purgati:
    for ripetizione in range(0,20):
        ricopia1=copy.deepcopy(network)
        random_disruption_relative_gc_dimension(ricopia1, 50, '065_random_disruption_relative_gc_dimension.csv')
