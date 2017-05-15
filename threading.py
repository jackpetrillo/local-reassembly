"""
Jack Petrillo and Rose Gold
Local Re-Assembly Final Project
CS 321

"""

import debruijn as db
import randomized_reads as rr
from random import randint
import copy as c


def thread(variant_reads, genome_db, genome_hash, gene_k, positions):
    """
    Takes as inputs a list of reads, the genome graph, the genome hash table, the
    kmer size in the graph, and the positions list of the reads.
    """
    ###GO THROUGH ALL READS
    for k in range(len(variant_reads)):
        read_kmers = db.kmerList(gene_k - 1, variant_reads[k])
        ###GO THROUGH KMERS OF SIZE GENE_K WITHIN READ
        for i in range(len(read_kmers)-1):
            cur_kmer = read_kmers[i] #current kmer in read
            if(cur_kmer in genome_db): #if already in the graph
                #check to see if it points to the next kmer in genome_db
                next_kmer = read_kmers[i+1]
                try:
                    index = genome_db[cur_kmer][0].index(next_kmer)
                except ValueError:
                    index = -1

                if(index > -1):
                    genome_db[cur_kmer][1][index] += 1
                else:
                    genome_db[cur_kmer][0].append(next_kmer)
                    genome_db[cur_kmer][1].append(1)


            else:
                genome_db[cur_kmer] = [[read_kmers[i+1]], [1]] #create a new node
                genome_hash[cur_kmer] = [positions[k] + i]

    return genome_db, genome_hash


def uniform_reads(mat, pat, ref, depth, read_len):
    """
    Takes as parameters maternal and paternal sequences, the reference gene, an average depth,
    and a read length.
    Returns num_kmers * depth amount of random reads, and their positions in the reference.
    """
    len_seq = len(ref)
    iterations = (len_seq - read_len + 1) * depth #to get average depth
    reads = []
    positions = []

    for i in range(iterations):
        mat_or_pat = randint(0, 2) #random stuff
        start_ind = randint(0, len_seq - read_len + 1)
        if mat_or_pat == 0:
            initial_read = mat[start_ind : start_ind + read_len]
            error_read = rr.rand_variant(initial_read, 30000) #adds errors
            reads.append(error_read)
            positions.append(start_ind)
        else:
            initial_read = pat[start_ind : start_ind + read_len]
            error_read = rr.rand_variant(initial_read, 30000) #adds errors
            reads.append(error_read)
            positions.append(start_ind)

    return reads, positions


def pruning(genome_db, threshold):
    """
    Takes as inputs a genome debruijn graph and an integer threshold. Prunes all edges
    of graphs that have weights under the threshold.
    """
    temp = c.deepcopy(genome_db) #not sure if needed -- make copy to ensure iteration works while deleting
    key_list = genome_db.keys()
    pruned = 0

    for index in range(len(key_list)):
        key = key_list[index]

        if(len(genome_db[key][1]) == 1): #if only one thing in list
            if(genome_db[key][1][0] < threshold): #and that thing has edge weight < threshold
                del temp[key] #delete it all
                pruned += 1

        elif(len(genome_db[key][1]) > 1): #if theres more...search
            for j in range(len(genome_db[key][1])):
                if(genome_db[key][1][j] < threshold): #remove bad stuff
                    temp[key][0].remove(genome_db[key][0][j])
                    temp[key][1].remove(genome_db[key][1][j])
                    pruned += 1


    return temp, pruned



def main():
    """
    Main method. Creates two sample sequences with SNVs. Gets reference material from debruijn.py.
    Creates random reads, threads them through the graph, and prunes the graph. Writes to files
    for debugging.
    """

    gene_k = 45 #This kmer size is huge. However, it was needed to avoid cycles in the DB graph, which we had a hard time getting around.
    prob = 50
    read_length = 50
    depth = 15 #average read depth
    prune_threshold = 10
    pruned = 0

    genome_db, kmer_positions, gene = db.main(gene_k) #gets reference stuff

    ###CREATES MATERNAL AND PATERNAL SNV SEQUENCES
    mat = gene[:-3] + 'Z' + gene[-2:]
    pat = gene[:-50] + "X" + gene[-49:]

    reads, positions = uniform_reads(mat, pat, gene, depth, read_length) #get list of reads, positions of reads.

    genome_db, genome_hash = thread(reads, genome_db, kmer_positions, gene_k, positions) #thread all reads

    genome_db, pruned = pruning(genome_db, prune_threshold) #prune

    #writes to files for debugging
    with open("genome_hash.txt", "w") as text_file:
        text_file.write(str(genome_hash))
    with open("genome_db.txt", "w") as text_file:
        text_file.write(str(genome_db))


    return genome_db, kmer_positions, gene, pruned


if __name__ == '__main__':
    main()
