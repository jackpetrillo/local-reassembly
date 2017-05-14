"""
Jack Petrillo and Rose Gold
Local Re-Assembly Final Project
CS 321

"""
import threading as th
import collections as c
import numpy as np

def path_scores(genome_db, genome_hash, ref_gene, start_kmer):
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



    #while start_node[0] in genome_db: #while the node is in the graph (Not the last node)
    # while (i < len(genome_db.keys())) and (start_node[0] in genome_db): #arbitrary? shouldnt be going through keys multiple times
    #
    #     node = start_node[0]
    #     denominator = sum(genome_db[node][1])
    #     print("The start node is at: " + str(genome_hash[node]))
    #
    #     for j in range(len(genome_db[node][0])):
    #         end_node = genome_db[node][0][j]
    #         print("It points to nodes at : " + str(genome_hash[end_node]))
    #
    #         path_prob = start_node[1]*(genome_db[node][1][j] / float(denominator))
    #
    #         to_search.append((end_node, path_prob))
    #
    #         if node in path_scores:
    #             path_scores[node][0].append(end_node)
    #             path_scores[node][1].append(path_prob)
    #         else:
    #             path_scores[node] = [[end_node],[path_prob]]
    #
    #     start_node = to_search.popleft()
    #
    #     i +=1



    #print(path_scores)

    # with open("likelihood_graph.txt", "w") as text_file:
    #     text_file.write(str(path_scores))






def main():
    genome_db, kmer_positions, gene = th.main()


    path_scores(genome_db, kmer_positions, gene, "CGTGACCCTCAGGTGATGCGCCAGGGCCGGCTGCCGTCGGGGAC")



if __name__ == '__main__':
    main()
