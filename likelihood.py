"""
Jack Petrillo and Rose Gold
Local Re-Assembly Final Project
CS 321

"""
import threading as th
import collections as c

def path_scores(genome_db, start_kmer):
    """
    Traverses all paths and makes a new graph with each paths likelihood scores.
    scores are computed as such:
    This score is calculated as the product of transition probabilities of the path edges,
    where the transition probability of an edge is computed as the number of reads supporting that
    edge divided by the sum of the support of all edges that share that same source vertex.
    """
    path_scores = {}

    start_node = (start_kmer, 1.0) #where we start...
    to_search = c.deque()

    while start_node[0] in genome_db: #while the node is in the graph (Not the last node)

        print("RUNNING")
        node = start_node[0]
        denominator = sum(genome_db[node][1])

        for j in range(len(genome_db[node][0])):
            end_node = genome_db[node][0][j]
            path_prob = start_node[1]*(genome_db[node][1][j] / float(denominator))
            to_search.append((end_node, path_prob))
            print("Start "  + node)
            print("END " + end_node)
            print("Score " + str(path_prob))

        start_node = to_search.popleft()


    print(path_scores)






def main():
    genome_db, kmer_positions, gene = th.main()
    if "CGTGACCCTCAGGT" in genome_db:
        print("THIS WORKED")


    path_scores(genome_db, "CGTGACCCTCAGGTGATGCGCCAGGGCCG")



if __name__ == '__main__':
    main()
