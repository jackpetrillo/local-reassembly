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

    node = start_kmer #where we start...
    to_search = c.deque()

    while node in genome_db: #while






def main():
    genome_db, kmer_positions, gene = th.main()


    path_scores(genome_db, "CGTGACCCTCAGGT")



if __name__ == '__main__':
    main()
