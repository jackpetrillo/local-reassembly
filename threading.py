"""
Jack Petrillo and Rose Gold
Local Re-Assembly Final Project
CS 321

"""

import debruijn as db
import randomized_reads as rr


def thread(variant_reads, genome_db, genome_hash, gene_k):
    """

    """

    ###GO THROUGH ALL READS
    for read in variant_reads:
        read_kmers = db.kmerList(gene_k - 1, read)
        ###GO THROUGH KMERS OF SIZE GENE_K WITHIN READ
        for i in range(len(read_kmers)-1):
            cur_kmer = read_kmers[i] #current kmer in read
            if(cur_kmer in genome_db): #if already in the graph
                #check to see if it points to the next kmer in genome_db
                for j in range(len(genome_db[cur_kmer][0])): #for every kmer this node points to
                    if(genome_db[cur_kmer][0][j] == read_kmers[i+1]): #if next in variant matches
                        genome_db[cur_kmer][1][j] += 1 #increment weight

            else:
                genome_db[cur_kmer] = [[read_kmers[i+1]], [1]] #create a new node




    return genome_db, genome_hash



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
    prob = 50
    read_length = 50

    genome_db, kmer_positions, gene = db.main(gene_k)

    variant1 = rr.rand_variant(gene, prob) #full sequence with SNVs w. probability prob
    variant1_reads = db.kmerList(read_length, variant1)


    genome_db, genome_hash = thread(variant1_reads, genome_db, kmer_positions, gene_k)

    variant2 = rr.rand_variant(gene, prob) #full sequence with SNVs w. probability prob
    variant2_reads = db.kmerList(read_length, variant2)


    genome_db, genome_hash = thread(variant2_reads, genome_db, kmer_positions, gene_k)



    print(genome_db)



if __name__ == '__main__':
    main()
