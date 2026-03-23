# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 20:27:59 2023

@author: lion
"""

import networkx, matplotlib.pyplot as plot

# Exercise 1
def kmers(sequence, k):
    n = len(sequence)
    kmers_list = []
    for i in range(n - k + 1):
        kmers_list += [sequence[i:i+k]]
    return kmers_list

# Exercise 2
def DeBruijnGraph(reads, k): # Steps in lecture 4 (slide 31)
    all_kmers = []
    for read in reads: # Step 1(b)
        all_kmers += kmers(read, k)

    graph = {'nodes':[], 'edges':[]} # Step 2

    for kmer in all_kmers: # Step 3
        prefix, suffix = kmer[:-1], kmer[1:]
        if prefix not in graph['nodes']: # (a)
            graph['nodes'] += [prefix]
        if suffix not in graph['nodes']: # (b)
            graph['nodes'] += [suffix]
        graph['edges'] += [[prefix, suffix]] # (c)

    return graph

def visualizeDBGraph(graph):
    dbGraph = networkx.DiGraph()
    
    dbGraph.add_nodes_from(graph['nodes']) #Add the nodes to the graph
    dbGraph.add_edges_from(graph['edges']) #Add the edges to the graph
    
    networkx.draw(dbGraph, with_labels=True, node_size=1000)
    plot.show()


        
#main()        
        
def build_degrees(graph, _in):
    assert _in == 1 or _in == 0, "_in must be 0 or 1"
    degrees = {}
    for node in graph['nodes']:
        degrees[node] = 0
        
    for node in graph['nodes']:
        for edge in graph['edges']:
            if node == edge[_in]: # When _in = 0, we are building the out-degrees, so we check that the node is a prefix by node == edge[_in]
                degrees[node] += 1
    return degrees

def buildAdjacencyList(graph):
    adjacencyList = {}
    for node in graph['nodes']:
        adjacencyList[node] = []

    for node in graph['nodes']:
        for edge in graph['edges']:
            if node == edge[0]: # We check that the node is a prefix by node == edge[0]
                adjacencyList[node] += [edge[1]] # We add its suffix to its adjacency list

    return adjacencyList

# Exercise 1
def GetStartNode(graph):
    _in = build_degrees(graph, _in = 1)
    out = build_degrees(graph, _in = 0) # _in = 0 means that we are building the out-degrees
    for node in graph['nodes']:
        if _in[node] == out[node] - 1:
            return node

def EulerianPath(graph): # Steps in lecture 5 (slide 7)        
        adjacencyList = buildAdjacencyList(graph)
        pathList = [] # Step 1
        startNode = GetStartNode(graph) # Step 2
        currentNode = startNode # Step 3

        def DFS(v):
            for u in adjacencyList[v]: # a
                adjacencyList[v].remove(u) # i
                DFS(u) # ii
            nonlocal pathList
            pathList += [v] # b
            
        DFS(currentNode) # Step 4
        return pathList

# Exercise 2
def AssembleGenome(pathList): # Steps in lecture 5 (slide 7)
    pathList = pathList[::-1] # Step 1
    genome = pathList[0] # Step 2
    for i in range(1,len(pathList)):
        genome += pathList[i][-1] # Step 3
    return genome
        
   

#main()  
        
        


k = None

# Exercise
def mergeChains(graph,k): # Lecture 6 (slide 17)
    _in = build_degrees(graph, _in = 1)
    out = build_degrees(graph, _in = 0)
    canMerge = True
    
    while canMerge:
        canMerge = False
        for edge in graph['edges']: # Condition 1
            A = edge[0] 
            B = edge[1]
            if out[A] == _in[B] == 1: # Condition 2
                canMerge = True
                #Merge
                graph['edges'].remove(edge)
                graph['nodes'].remove(A); graph['nodes'].remove(B)
                newNode = A + B[k-1:]
                graph['nodes'].append(newNode) # Update nodes
                
                for e in graph['edges']: # Update edges
                    if e[0] == B:
                        e[0] = newNode
                    if e[1] == A:
                        e[1] = newNode
                
                _in[newNode] = _in[A]; out[newNode] = out[B] # Update in and out degrees
                _in.pop(A); _in.pop(B)
                out.pop(A); out.pop(B)
                
    return graph

def main():
    global k
    k = 4
    graph = DeBruijnGraph(['TTACGTT','CCGTTA','GTTAC','GTTCGA','CGTTC'],5)
    # visualizeDBGraph(graph)
    newGraph = mergeChains(graph,k)
    # visualizeDBGraph(newGraph)
    print(newGraph)
    # Example in lecture 6 (slide 18)
    k = 3
    # graph2 = {'nodes':['GAC', 'ACC', 'CCG', 'CGT', 'GTA', 'TAA', 'AAT', 'TAG', 'ACT', 'CTG', 'TGT'],
    #           'edges':[['GAC', 'ACC'], ['GAC', 'ACT'], ['ACC', 'CCG'], ['CCG', 'CGT'], ['CGT', 'GTA'], ['GTA', 'TAA'], ['TAA', 'AAT'],
    #                    ['GTA', 'TAG'], ['ACT', 'CTG'], ['CTG','TGT'], ['TGT', 'GTA']]}
    # visualizeDBGraph(graph2)
    # newGraph = mergeChains(graph2,k)
    # visualizeDBGraph(newGraph)
    print(newGraph)
main()     
