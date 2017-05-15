"""
Jack Petrillo and Rose Gold
Local Re-Assembly Final Project
CS 321

"""
import threading as th
import collections as c
import numpy as np

def path_score(genome_db, genome_hash, ref_gene, start_kmer):
    """
    Traverses all paths and makes a new graph with each paths likelihood scores.
    scores are computed as such:
    This score is calculated as the product of transition probabilities of the path edges,
    where the transition probability of an edge is computed as the number of reads supporting that
    edge divided by the sum of the support of all edges that share that same source vertex.
    """
    keys = genome_db.keys()
    kmer_size = len(keys[0])
    num_kmers = len(ref_gene) - kmer_size + 1

    path_kmers = [[] for i in range(num_kmers)]
    path_scores = [[] for i in range(num_kmers)]


    start_node = (start_kmer, 1.0) #where we start...
    to_search = c.deque()


    while(start_node[0] in genome_db):

        node = start_node[0]
        node_position = genome_hash[node][0] #picking arbitrary first location
        node_prob = start_node[1]
        denominator = sum(genome_db[node][1])


        if node in path_kmers[node_position]:
            index = path_kmers[node_position].index(node) #get the index
            path_scores[node_position][index] += node_prob #at that index, update scores with its probability
        else:
            path_kmers[node_position].append(node)
            path_scores[node_position].append(node_prob)

        node_index = path_kmers[node_position].index(node) #this will always be a find.

        for j in range(len(genome_db[node][0])): #for everything it points to
            end_node = genome_db[node][0][j] #get the child
            path_prob = path_scores[node_position][node_index] * (genome_db[node][1][j] / float(denominator))

            to_search.append((end_node, path_prob))

        start_node = to_search.popleft()
        
    return path_kmers, path_scores


def variant_regions(path_kmers, path_scores):
    """
    Returns two arrays, where each indices corresponds to a start index and end index
    for potential variant sites.
    """
    var_starts = []
    var_ends = []

    run = False

    for i in range(len(path_kmers)):
        if(run == False and len(path_kmers[i]) > 1):
            run = True
            var_starts.append(i)

        if(run == True and len(path_kmers[i]) == 1):
            var_ends.append(i)
            run = False

        if(run == True and i == len(path_kmers) - 1):
            var_ends.append(i)


    return var_starts, var_ends


def main():
    """
    Runs main method of threading.py to get the threaded graph, hash table, gene,
    and # of pruned edges.
    Passes to path_score() which calculates path proabibility scores.
    Passes to variant_regions() which traces the graph and returns arrays with
    locations of possible variants.
    Prints results.
    """
    genome_db, kmer_positions, gene, pruned = th.main()

    #hard coded for 45
    path_kmers, path_scores = path_score(genome_db, kmer_positions, gene, "CGTGACCCTCAGGTGATGCGCCAGGGCCGGCTGCCGTCGGGGAC")
    #gets arrays for variant stars/ends
    var_starts, var_ends = variant_regions(path_kmers, path_scores)
    #prints pruned edges
    print("Pruned edges: " + str(pruned))
    #Prints variants
    for j in range(len(var_ends)):
        print(str(j) + ": Potential Variant Start: " + str(var_starts[j]) + "\t Variant End: " +  str(var_ends[j]))

    #writes to files for debugging
    with open("haplo_likelihood.txt", "w") as text_file:
        for j in range(len(var_ends)):
            text_file.write(str(var_starts[j]) + '\n')
            text_file.write(str(var_ends[j]) + '\n')


if __name__ == '__main__':
    main()
