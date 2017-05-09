"""
Jack Petrillo
Rosalind 3/16 - 'Construct the De Bruijn Graph of a String'
CS 321

"""

import sys
import kmer


def dict_create(kmerList):
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


if __name__ == "__main__":

    """
    The main method reads an integer k and a string from an input file. Passes the list to
    the dict_create() function which returns a dictionary of an adjacency list representation
    of the overlap graph.
    """

    f = open(sys.argv[1], "r")

    strk = f.readline().strip()
    k = (int)(strk) ##typecast to integer
    word = f.readline().strip()

    klist = kmer.kmerList(k-1, word) #create kmer list (with k - 1)
    #print(klist)

    final_dict = dict_create(klist) #adjacency list dictionary

    #print(final_dict)
    for key in final_dict: #printing all node -> outgoing nodes
        outgoing = ""

        for i in range(len(final_dict[key])):
            outgoing = outgoing + final_dict[key][i] + ","

        outgoing = outgoing[:-1]
        
        print(key + " -> " + outgoing)
