"""
Jack Petrillo
Rosalind 3/16 - 'Generate the k-mer Composition of a String'
CS 321

"""

import sys

def kmerList(k, word):
    """
    Takes as inputs an integer k and a String word. Returs a list of all the k-mers of word.
    """

    klist = []
    for index in range(len(word) - k + 1):
        klist.append(word[index:index+k])

    return klist

if __name__ == "__main__":
    """
    The main method reads an integer k and a string from an input file.
    Prints all of the k-mers of the string.
    """

    f = open(sys.argv[1], "r")

    strk = f.readline().strip()
    k = (int)(strk) ##typecast to integer
    word = f.readline().strip()

    klist = kmerList(k, word)

    for element in klist: #print all elements of list
        print(element)
