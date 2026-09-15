#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import networkx as nx
import copy



    
def purge_smaller_components(network, minimum_number_of_nodes_permitted): #minimum must be int and not lower than 2
    print(network)
    sorted_components = sorted(nx.connected_components(network), key=len, reverse=True)
    print(len(sorted_components))
    if len(sorted_components)==1:
        purged_network=network.subgraph(network.nodes)
    else:
        for i in range(0, len(sorted_components)):
            print(len(sorted_components[i]))
            if len(sorted_components[i])<minimum_number_of_nodes_permitted:
                purged_network=network.subgraph(set().union(*sorted_components[0:i]))
                break
    return purged_network

def sort_nodes_by_degree(network):
    sorted_nodes=sorted(nx.degree(network), key = lambda x: x[1], reverse=True) #ma ce devo mette .items o no?
    return sorted_nodes


def sort_nodes_by_katz_centrality(network):
    sorted_nodes=sorted(nx.katz_centrality_numpy(network, alpha=0.06).items(), key = lambda x: x[1], reverse=True)
    return sorted_nodes

def sort_nodes_by_betweenness_centrality(network):
    sorted_nodes=sorted(nx.betweenness_centrality(network).items(), key = lambda x: x[1], reverse=True)
    return sorted_nodes

def remove_nodes_by_degree_and_calculate_absolute_efficiency(unfrozen_network, number_of_iterations, filename):
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
                f.write(sort_nodes_by_degree(unfrozen_network)[0][0]+':')
                unfrozen_network.remove_node(sort_nodes_by_degree(unfrozen_network)[0][0])
        
def remove_nodes_by_degree_and_calculate_relative_efficiency(unfrozen_network, number_of_iterations, filename):
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
                f.write(sort_nodes_by_degree(unfrozen_network)[0][0]+':')
                unfrozen_network.remove_node(sort_nodes_by_degree(unfrozen_network)[0][0])

def remove_nodes_by_degree_and_calculate_absolute_gc_dimension(unfrozen_network, number_of_iterations, filename):
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
                f.write(sort_nodes_by_degree(unfrozen_network)[0][0]+':')
            unfrozen_network.remove_node(sort_nodes_by_degree(unfrozen_network)[0][0])

def remove_nodes_by_degree_and_calculate_relative_gc_dimension(unfrozen_network, number_of_iterations, filename):
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
                f.write(sort_nodes_by_degree(unfrozen_network)[0][0]+':')
                unfrozen_network.remove_node(sort_nodes_by_degree(unfrozen_network)[0][0])
    
    

def remove_nodes_by_katz_centrality_and_calculate_absolute_efficiency(unfrozen_network, number_of_iterations, filename):
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
                f.write(sort_nodes_by_katz_centrality(unfrozen_network)[0][0]+':')
                unfrozen_network.remove_node(sort_nodes_by_katz_centrality(unfrozen_network)[0][0])
        
def remove_nodes_by_katz_centrality_and_calculate_relative_efficiency(unfrozen_network, number_of_iterations, filename):
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
                f.write(sort_nodes_by_katz_centrality(unfrozen_network)[0][0]+':')
                unfrozen_network.remove_node(sort_nodes_by_katz_centrality(unfrozen_network)[0][0])

def remove_nodes_by_katz_centrality_and_calculate_absolute_gc_dimension(unfrozen_network, number_of_iterations, filename):
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
                f.write(sort_nodes_by_katz_centrality(unfrozen_network)[0][0]+':')
            unfrozen_network.remove_node(sort_nodes_by_katz_centrality(unfrozen_network)[0][0])

def remove_nodes_by_katz_centrality_and_calculate_relative_gc_dimension(unfrozen_network, number_of_iterations, filename):
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
                f.write(sort_nodes_by_katz_centrality(unfrozen_network)[0][0]+':')
                unfrozen_network.remove_node(sort_nodes_by_katz_centrality(unfrozen_network)[0][0])
    



    
    
    
def remove_nodes_by_betweenness_and_calculate_absolute_efficiency(unfrozen_network, number_of_iterations, filename):
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
                f.write(sort_nodes_by_betweenness_centrality(unfrozen_network)[0][0]+':')
                unfrozen_network.remove_node(sort_nodes_by_betweenness_centrality(unfrozen_network)[0][0])
        
def remove_nodes_by_betweenness_and_calculate_relative_efficiency(unfrozen_network, number_of_iterations, filename):
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
                f.write(sort_nodes_by_betweenness_centrality(unfrozen_network)[0][0]+':')
                unfrozen_network.remove_node(sort_nodes_by_betweenness_centrality(unfrozen_network)[0][0])

def remove_nodes_by_betweenness_and_calculate_absolute_gc_dimension(unfrozen_network, number_of_iterations, filename):
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
                f.write(sort_nodes_by_betweenness_centrality(unfrozen_network)[0][0]+':')
                unfrozen_network.remove_node(sort_nodes_by_betweenness_centrality(unfrozen_network)[0][0])

def remove_nodes_by_betweenness_and_calculate_relative_gc_dimension(unfrozen_network, number_of_iterations, filename):
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
                f.write(sort_nodes_by_betweenness_centrality(unfrozen_network)[0][0]+':')
                unfrozen_network.remove_node(sort_nodes_by_betweenness_centrality(unfrozen_network)[0][0])


def write_clustering_and_efficiency(networks, minimum_number_of_nodes_permitted):
    for i in range(1, minimum_number_of_nodes_permitted+1):
        with open('clustering_and_efficiency.txt', 'a') as f:
            if i==1:
                f.write("nessun nodo escluso"+'\n')
                for network in networks:
                    f.write(network.name+' clustering coefficient = '+str(nx.average_clustering(network))+'\n')
                    f.write(network.name+' global efficiency = '+str(nx.global_efficiency(network))+'\n')
            else:
                f.write("esclusi componenti sotto i "+str(i)+" nodi\n")
                for network in networks:
                    purged_network=purge_smaller_components(network, i)
                    f.write(network.name+' clustering coefficient = '+str(nx.average_clustering(purged_network))+'\n')
                    f.write(network.name+' global efficiency = '+str(nx.global_efficiency(purged_network))+'\n')
                    

def write_katz_centrality(networks, filename):
   with open(filename, 'a') as f:
        for network in networks:
            f.write(network.name+'\n') 
            sorted_nodes=sort_nodes_by_katz_centrality(network)
            for node in sorted_nodes:
                f.write(node[0] +' = '+str(node[1])+'\n')
            f.write('\n\n\n')
    
def write_betweenness_centrality(networks, filename):
   with open(filename, 'a') as f:
        for network in networks:
            f.write(network.name+'\n') 
            sorted_nodes=sort_nodes_by_betweenness_centrality(network)
            for node in sorted_nodes:
                f.write(node[0] +' = '+str(node[1])+'\n')
            f.write('\n\n\n')
    
    

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

networks_065=[HC, MS, MC, SS, SC]
networks_purgati_065=[]


for network in networks_065:
    print(network)
    nx.relabel_nodes(network, nx.get_node_attributes(network, 'shared name'), copy=False)
    network_purgato=purge_smaller_components(network, 5)
    network_purgato_scongelato=nx.Graph(network_purgato)
    networks_purgati_065.append(network_purgato_scongelato)


for network in networks_purgati_065:
    ricopia1=copy.deepcopy(network)
    remove_nodes_by_betweenness_and_calculate_relative_gc_dimension(ricopia1, 30, '065_relative_gc_dimension_betweenness.csv')

for network in networks_purgati_065:
    ricopia1=copy.deepcopy(network)
    remove_nodes_by_katz_centrality_and_calculate_relative_gc_dimension(ricopia1, 30, '065_relative_gc_dimension_katz_centrality.csv')

for network in networks_purgati_065:
    ricopia1=copy.deepcopy(network)
    remove_nodes_by_degree_and_calculate_relative_gc_dimension(ricopia1, 30, '065_relative_gc_dimension_degree.csv')


