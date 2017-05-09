"""
Jack Petrillo and Rose Gold
Local Re-Assembly Final Project
CS 321

"""

import sys

def get_gene(filename):

    with open(filename, 'r') as file_object:
        file_object.readline()
        input_str = file_object.read().strip()
    return input_str


def create_Debruijn(kmerList):
    """
    Takes as an input a list of kmers. Returns a dictionary of the adjacency list
    representation of the de bruijn graph.
    """
    dictionary = {}

    for i in range(len(kmerList)-1): #iterate over all kmers (besides last)
        if(kmerList[i] in dictionary):
            templist = dictionary[kmerList[i]] #take list
            templist.append(kmerList[i+1]) #add new outgoing edge
            dictionary[kmerList[i]] = templist #update list
        else:
            dictionary[kmerList[i]] = [kmerList[i+1]] #create 1 element list

    return dictionary


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
    """

    for key in graph: #printing all node -> outgoing nodes
        outgoing = ""

        for i in range(len(graph[key])):
            outgoing = outgoing + graph[key][i] + ","

        outgoing = outgoing[:-1]

        print(key + " -> " + outgoing)


def main():
    """
    The main method reads a string from an input file. Creates a kmer list and passes the list to
    the dict_create() function which returns a dictionary of an adjacency list representation
    of the overlap de bruijn graph.
    """

    gene = get_gene("OPN1LW.txt")

    k = 10 #for now set to constant

    klist = kmer.kmerList(k-1, gene) #create kmer list (with k - 1)

    final_dict = create_Debruijn(klist) #adjacency list dictionary

    print_Debruijn(final_dict)


if __name__ == "__main__":
    main()
