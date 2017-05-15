"""
Jack Petrillo and Rose Gold
Local Re-Assembly Final Project
CS 321

"""

import sys

def get_gene(filename):
    """
    Takes as an input a file name to read.
    Ignores the first line of the file (has identifying info in applicable genome files).
    Returns a string of gene devoid of new line chars.
    """
    with open(filename, 'r') as file_object:
        file_object.readline()
        input_str = file_object.read().strip()
        input_str = input_str.replace("\n", "")
    return input_str


def create_Debruijn(kmerList):
    """
    Takes as an input a list of kmers. Returns a dictionary of the adjacency list
    representation of the de bruijn graph.
    """
    dictionary = {}
    kmer_positions = {} #hash table of unique kmers

    for i in range(len(kmerList)-1): #iterate over all kmers (besides last)
        if(kmerList[i] in dictionary):
            templist = dictionary[kmerList[i]] #take list
            templist[0].append(kmerList[i+1]) #add new outgoing edge
            templist[1].append(0)
            dictionary[kmerList[i]] = templist #update list

            kmer_positions[kmerList[i]].append(i) #update locations
        else:
            dictionary[kmerList[i]] = [[kmerList[i+1]], [0]] #create 1 element list
            kmer_positions[kmerList[i]] = [i] #create hash table location

    return dictionary, kmer_positions


def kmerList(k, word):
    """
    Takes as inputs an integer k and a String word. Returs a list of all the k-mers of word.
    """

    klist = []
    for index in range(len(word) - k + 1):
        klist.append(word[index:index+k])

    return klist


def print_Debruijn(graph):
    """
    Passed an adjacency list representation of a de bruijn graph. Prints the graph.
    Used for debugging.
    """

    for key in graph: #printing all node -> outgoing nodes
        outgoing = ""

        for i in range(len(graph[key][0])):
            outgoing = outgoing + graph[key][0][i] + " " + str(graph[key][1][i]) + ","

        outgoing = outgoing[:-1]

        print(key + " -> " + outgoing)


def main(k):
    """
    Takes an integer k as a parameters.
    The main method reads a string from an input file.  Creates a kmer list and passes the list to
    the dict_create() function which returns a dictionary of an adjacency list representation
    of the overlap de bruijn graph, as well as the kmer position hash table, and the gene.
    """

    gene = get_gene("OPN1LW.txt")

    klist = kmerList(k-1, gene) #create kmer list (with k - 1)

    final_dict, kmer_positions = create_Debruijn(klist) #adjacency list dictionary

    return (final_dict, kmer_positions, gene)



if __name__ == "__main__":
    main()
