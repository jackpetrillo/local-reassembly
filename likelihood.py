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


    for index in range(len(path_kmers)):
        print("Position: " + str(index))
        for j in range(len(path_kmers[index])):
            print("kmer: " + path_kmers[index][j] + " with probability " + str(path_scores[index][j]))

    return path_kmers, path_scores


def variant_regions(path_kmers, path_scores):

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
    genome_db, kmer_positions, gene = th.main()


    path_kmers, path_scores = path_score(genome_db, kmer_positions, gene, "CGTGACCCTCAGGTGATGCGCCAGGGCCGGCTGCCGTCGGGGAC")

    var_starts, var_ends = variant_regions (path_kmers, path_scores)

    for j in range(len(var_ends)):
        print("variant start: " + str(var_starts[j]))
        print("variant end: " + str(var_ends[j]))



if __name__ == '__main__':
    main()
