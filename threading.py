"""
Jack Petrillo and Rose Gold
Local Re-Assembly Final Project
CS 321

"""

import debruijn as db



def thread(variant, genome_db, genome_hash):
    """

    """





def main():
    """
    Main method. Has a sample variant and k. Reads genome from file, creates, a kmerlist,
    creates the genome graph and kmer hash. Threads the sample variant into the db graph.
    """
    # change this stuff later
    ##50 chars...ending 3 chars from end of genome. Changing last 3 from AAC to TTT
    sample_variant = "CCCTCCTTCTCCATCCCTGTAAAATAAATGTAATTTATCTTTGCCAATTT"
    #constant k
    k = 25

    gene = db.get_gene("OPN1LW.txt")
    genome_kmerlist = db.kmerList(k, gene)
    genome_db, kmer_positions = db.create_Debruijn(genome_kmerlist)

    thread(sample_variant, genome_db, kmer_positions)




if __name__ == '__main__':
    main()
