from collections import defaultdict
# import networkx, matplotlib.pyplot as plot



def kmers(seq, k):
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]


def DeBruijnGraph(reads, k):
    dicti = {}
    E = []
    for read in reads:
        KMers = kmers(read, k)
        edges = [kmers(mer, k - 1) for mer in KMers]
        for edge in edges:
            if edge[0] not in dicti.keys(): dicti[edge[0]] = []
            if edge[1] not in dicti.keys(): dicti[edge[1]] = []
            dicti[edge[0]].append(edge[1])
            E.append(edge)
    V = list(dicti.keys())
    graph = {'nodes': V, 'edges': E}
    print(graph)
    return (graph, dicti)



def generate_eulerian_path(dna_seq,k):
    graph,de_bruijn_graph = DeBruijnGraph(dna_seq,k)
    path = []
    stack = [list(de_bruijn_graph.keys())[len(list(de_bruijn_graph.keys()))-len(dna_seq)]]
    while stack:
        u = stack[-1]
        if de_bruijn_graph[u]:
            stack.append(de_bruijn_graph[u].pop())
        else:
            path.append(stack.pop())
    arr=[]
    arr.extend(path[::-1])
    return arr



dna_seq = ['TTACGTT', 'CCGTTA', 'GTTAC', 'GTTCGA', 'CGTTC']
graph,dicti=DeBruijnGraph(dna_seq,5)
eulerian_path = generate_eulerian_path(dna_seq,5)
print("Eulerian Path: ", eulerian_path)
