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
        
from Lecture_4 import DeBruijnGraph

def main():
    graph = DeBruijnGraph(['TTACGTT','CCGTTA','GTTAC','GTTCGA','CGTTC'],5)
    path = EulerianPath(graph)
    print(path)
    print(AssembleGenome(path))

#main()  
        
        
