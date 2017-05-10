"""
Jack Petrillo and Rose Gold
Local Re-Assembly Final Project
CS 321
Generate Randomized Reads from a Sequence
"""

from random import randint

def get_gene(filename):

    with open(filename, 'r') as file_object:
        file_object.readline()
        input_str = file_object.read().strip()
        input_str = input_str.replace("\n", "")
    return input_str
