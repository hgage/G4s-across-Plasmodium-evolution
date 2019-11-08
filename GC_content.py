"""This script calculates the size and GC content of whole genomes, as well as the GC content of var genes"""

from Bio import SeqIO

def calculateGC(seq):
    """Calculate the GC content of a DNA sequence

        :param seq- a string representation of a DNA sequence
        :return a floating point decimal of the % GC content in the sequence
    """
    sequence_length = 0  # the total number of nucleotides in the sequence
    GC = 0  # the number of G and C nucleotides in the sequence
    for base in str(seq):
        if base.upper() == "G" or base.upper() == "C":
            GC += 1
        sequence_length += 1
    GC_content = (GC / sequence_length) * 100
    return GC_content

# store the paths to the files containing annotated transcripts for each Laveranian species
Pfalciparum_path = "/Users/HunterGage/Desktop/Gene Annotations/PlasmoDB-44_Pfalciparum3d7_AnnotatedCDSs.txt"
Padleri_path = "/Users/HunterGage/Desktop/Gene Annotations/PlasmoDB-44_PadleriG01_AnnotatedCDSs.txt"
Pgaboni_path = "/Users/HunterGage/Desktop/Gene Annotations/PlasmoDB-44_PgaboniG01_AnnotatedCDSs.txt"
Pbillcollinsi_path = "/Users/HunterGage/Desktop/Gene Annotations/PlasmoDB-44_PbillcollinsiG01_AnnotatedCDSs.txt"
Pblacklocki_path = "/Users/HunterGage/Desktop/Gene Annotations/PlasmoDB-44_PblacklockiG01_AnnotatedCDSs.txt"
Ppraefalciparum_path = "/Users/HunterGage/Desktop/Gene Annotations/PlasmoDB-44_PpraefalciparumG01_AnnotatedCDSs.txt"
Preichenowi_path = "/Users/HunterGage/Desktop/Gene Annotations/PlasmoDB-44_PreichenowiG01_AnnotatedCDSs.txt"

species_list = [Pfalciparum_path, Padleri_path, Pgaboni_path, Pbillcollinsi_path, Pblacklocki_path, Ppraefalciparum_path,
           Preichenowi_path]

# populate dictionary, where key= species and value= a list of SeqIO objects with information about coding sequences
species_annotations = {}
for species in species_list:
    seqs_info = []  # list of SeqIO objects with information about coding sequences for a species
    # represent each sequence in the file of annotated transcripts as a SeqIO object, and add to list of SeqIO objects
    for seq_record in SeqIO.parse(species, "fasta"):
        seqs_info.append(seq_record)
    species_annotations.update({species.split("_")[1]: seqs_info})

def cleandescriptions(seqs_info):
    """split descriptions from SeqIO objects into a list and remove whitespace from the entries
            description[0]: gene ID
            description[1]: organism
            description[2]: gene annotation
            description[3]: location
            description[4]: transcript length

    :param seqs_info: a list of SeqIO objects containing information about coding sequences for a species
    :return: list of SeqIO objects, now with cleaned descriptions so the gene annotations can be easily accessed
    """

    for seq in seqs_info:
        # split description to list, with indices for geneID, organism, gene annotation, location, and transcript length
        seq.description = seq.description.split("|")
        seq.description = [x.strip(' ') for x in seq.description]  # strip whitespace
    return seqs_info

# clean the descriptions in the SeqIO objects for all of the species, so gene annotations can be easily accessed
for species in species_annotations:
    cleandescriptions(species_annotations[species])

# populate dictionary, where key= species and value= a list of var gene sequences for the species
species_var_genes = {}
for species in species_annotations:
    var_genes = []  # list of var gene sequences for the species
    for gene in species_annotations[species]:
        if "erythrocyte membrane protein" in gene.description[2].lower():
            var_genes.append(gene.seq)
    species_var_genes.update({species: var_genes})

# populate dictionary, where key= species and value = % GC content of var genes in the species
summary_data = {}
for species in species_var_genes:
    long_var_gene_seq = ""  # concatenation of all var gene sequences for the species
    for var_gene_sequence in species_var_genes[species]:
        long_var_gene_seq += var_gene_sequence
    summary_data.update({species: calculateGC(long_var_gene_seq)})

# print data for GC content of var genes for the Laveranian species
print("var gene GC contents: ")
for entry in summary_data:
    print(entry, summary_data[entry])

# store the paths to the files containing whole genome sequences for each species
Pfalciparum_genome_path = "/Users/HunterGage/Desktop/Whole Genome Sequences/PlasmoDB-44_Pfalciparum3d7_Genome.txt"
Padleri_genome_path = "/Users/HunterGage/Desktop/Whole Genome Sequences/PlasmoDB-44_PadleriG01_Genome.txt"
Pgaboni_genome_path = "/Users/HunterGage/Desktop/Whole Genome Sequences/PlasmoDB-44_PgaboniG01_Genome.txt"
Pblacklocki_genome_path = "/Users/HunterGage/Desktop/Whole Genome Sequences/PlasmoDB-44_PblacklockiG01_Genome.txt"
Pbillcollinsi_genome_path = "/Users/HunterGage/Desktop/Whole Genome Sequences/PlasmoDB-44_PbillcollinsiG01_Genome.txt"
Preichenowi_genome_path = "/Users/HunterGage/Desktop/Whole Genome Sequences/PlasmoDB-44_PreichenowiG01_Genome.txt"
Ppraefalciparum_genome_path = "/Users/HunterGage/Desktop/Whole Genome Sequences/PlasmoDB-44_PpraefalciparumG01_Genome.txt"
Prelictum_genome_path = "/Users/HunterGage/Desktop/Whole Genome Sequences/PlasmoDB-44_PrelictumSGS1-like_Genome.txt"
Pgallinaceum_genome_path = "/Users/HunterGage/Desktop/Whole Genome Sequences/PlasmoDB-44_Pgallinaceum8A_Genome.txt"
Pberghei_genome_path = "/Users/HunterGage/Desktop/Whole Genome Sequences/PlasmoDB-44_PbergheiANKA_Genome.txt"
Pchabaudi_genome_path = "/Users/HunterGage/Desktop/Whole Genome Sequences/PlasmoDB-45_Pchabaudichabaudi_Genome.txt"

species_genome_list = [Pfalciparum_genome_path, Padleri_genome_path, Pgaboni_genome_path, Pbillcollinsi_genome_path,
                       Pblacklocki_genome_path, Ppraefalciparum_genome_path, Preichenowi_genome_path, Prelictum_genome_path,
                       Pgallinaceum_genome_path, Pberghei_genome_path, Pchabaudi_genome_path]

# populate dictionary, where key= species and value= GC content of the whole genome for the species
summary_data_genomeGC = {}
for species in species_genome_list:
    seqs_info = []  # list of SeqIO objects containing information about genomic sequences for the species
    long_sequence = ""  # concatenation of all genomic sequences for the species
    genome_size = 0  # number of nucleotides in the whole genome
    for seq_record in SeqIO.parse(species, "fasta"):
        seqs_info.append(seq_record.seq)
    for seq in seqs_info:
        genome_size += len(str(seq))
    # print data for the sizes of the genomes
    print(species.split("_")[1], "whole genome size:", genome_size)
    for seq in seqs_info:
        long_sequence += str(seq)
    summary_data_genomeGC.update({species.split("_")[1]: calculateGC(long_sequence)})

# print data for GC content of whole genome for all species
print("Whole genome GC contents:")
for entry in summary_data_genomeGC:
    print(entry, summary_data_genomeGC[entry])







