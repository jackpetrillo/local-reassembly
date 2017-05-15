"""
Jack Petrillo and Rose Gold
Local Re-Assembly Final Project
CS 321
"""
from graphviz import Digraph
from string import ascii_uppercase

# from graphviz import render
# #import graphviz
#
# dot = Digraph(comment = 'The Round Table')
# dot.node("A", "King Arthur")
# dot.node("B", "Sir Bedevere the Wise")
# dot.node("L", "Sir Lancelot the Brave")
# dot.edges(["AB", "AL"])
# dot.edge("B", "L", constraint = "false")
#
# print(dot)
#
# dot.render('graphviz_trial.gv', view=True)

def get_likelihood_pos(filename):
    """
    """
    positions = []
    with open(filename, 'r') as file_object:
        for line in file_object:
            line = line.strip()
            positions.append(line)

    return positions


def output_digraph(positions, len_gene):
    """
    """
    node[style=filled, color=cornflowerblue, fontcolor=white, fontsize=10,
        fontname="Helvetica"]
    edge[arrowhead=vee, arrowtail=inv, arrowsize=.7, color=maroon, fontsize=10,
        fontcolor=navy]
    alph_ind = 0
    cur_char = ascii_uppercase[alph_ind]

    dot = Digraph(comment = "OPN1LW Haplotypes")
    dot.node(cur_char, str(alph_ind))

    for index in positions:
        next_char = ascii_uppercase[alph_ind + 1]
        dot.node(next_char, str(index))
        dot.edge(cur_char, next_char, constraint = "false")
        if alph_ind % 2 != 0:
            dot.edge(cur_char, next_char, constraint = "false")
        cur_char = ascii_uppercase[alph_ind + 1]
        alph_ind += 1

    dot.node("Z", str(len_gene - 1))
    dot.edge(cur_char, "Z", constraint = "false")
    return dot

def main():
    """
    """
    filename = "haplo_likelihood.txt"
    len_gene = 14808
    positions = get_likelihood_pos(filename)
    dot = output_digraph(positions, len_gene)
    print(dot)

if __name__ == '__main__':
    main()
