"""
Jack Petrillo and Rose Gold
Local Re-Assembly Final Project
CS 321
Generate Randomized Reads from a Sequence
"""

def kmerList(k, word):
    """
    Takes as inputs an integer k and a String word. Returs a list of all the k-mers of word.
    """

    klist = []
    for index in range(len(word) - k + 1):
        klist.append(word[index:index+k])

    return klist
