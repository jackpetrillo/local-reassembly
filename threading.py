"""
Jack Petrillo and Rose Gold
Local Re-Assembly Final Project
CS 321

"""

import debruijn as db
import randomized_reads as rr


def thread(read, genome_db, genome_hash, gene_k):
    """

    """
    read_kmers = db.kmerList(gene_k, read)

    cur_kmer = read_kmers[0]
    start_position = -1

    if cur_kmer in genome_hash:
        start_position = genome_hash[cur_kmer]
    else:
        genome_hash[cur_kmer] = 0
        genome_db[cur_kmer] = [[read_kmers[1]],[1]]





def get_Variant(genome, k):
    """

    """
    var_seq = rr.rand_variant(genome, 10) #10 is probability
    read_list = db.kmerList(k, var_seq)

    return read_list[-1] #returns last element of read_list


def main():
    """
    Main method. Has a sample variant and k. Reads genome from file, creates, a kmerlist,
    creates the genome graph and kmer hash. Threads the sample variant into the db graph.
    """
    # gene_k = 25
    # gene = db.get_gene("OPN1LW.txt")
    # genome_kmerlist = db.kmerList(gene_k, gene)
    # genome_db, kmer_positions = db.create_Debruijn(genome_kmerlist)

    gene_k = 25

    genome_db, kmer_positions, gene = db.main(gene_k)

    read_k = 50
    sample_variant = get_Variant(gene, read_k)

    #print(gene)
    #print(sample_variant)

    thread(sample_variant, genome_db, kmer_positions, gene_k)




if __name__ == '__main__':
    main()
