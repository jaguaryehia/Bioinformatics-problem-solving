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

def main():
    graph = DeBruijnGraph(['TTACGTT','CCGTTA','GTTAC','GTTCGA','CGTTC'],5)
    visualizeDBGraph(graph)
        
#main()        
        
