"""This script tests for a statistical co-distribution of vars/pirs with PQSs. It calculates the mean and median
distance of each var/pir gene from its nearest PQS in the actual data set, and compares these values to a data set
consisting of 1,000,000 simulated var/pir genes with sizes normally distributed around the mean var/pir gene size.
"""

# import some useful libraries
from Bio import SeqIO
from collections import defaultdict
import csv
from matplotlib import pyplot as plt
import numpy
import random
from scipy import stats
import statistics

# list of sequence record objects for each gene annotation in the fasta file
seqs_info = []
for seq_record in SeqIO.parse("/Users/HunterGage/Desktop/Gene Annotations/PlasmoDB-44_Pfalciparum3d7_AnnotatedCDSs.txt",
                              "fasta"):
    seqs_info.append(seq_record)

def cleandescriptions():
    """Split descriptions from SeqIO objects into a list and remove whitespace from the entries, so the information can
    be easily accessed.
            description[0]: gene ID
            description[1]: organism
            description[2]: gene annotation
            description[3]: location
            description[4]: transcript length
    """
    for seq in seqs_info:
        seq.description = seq.description.split("|")
        seq.description = [x.strip(' ') for x in seq.description]
    return seqs_info

cleandescriptions()

def cleanlocations():
    """Split the location in the SeqRecord into a list so that the values for "start", "end", and "strand" can be easily
    accessed.
        to access start location of coding sequence: seq.description[3][1][0]
        to access end location of coding sequence: seq.description[3][1][1]
        to access strand: seq.description[3][2][:seq.description[3][2].find("\)")]
    """
    for seq in seqs_info:
        seq.description[3] = seq.description[3].replace("(", ":")
        seq.description[3] = seq.description[3].split(":")
        seq.description[3][1] = seq.description[3][1].split("-")
        seq.description[3][1] = list(map(int, seq.description[3][1]))
    return seqs_info

cleanlocations()

PQS_dict = defaultdict(list)  # stores data for PQSs on the sense strand
PQS_dict_RC = defaultdict(list)  # stores data for PQSs on the reverse-complement strand
def cleanPQSfile():
    """Read from csv file containing PQSs and their locations. Store PQS location, chromosome, and length in PQS_dict
    and PQS_dict_RC.
    Dictionary structure: key = PQS start location
                          value = a list, where index 0 = chromosome number, and index 1 = PQS length
    """
    with open("/Users/HunterGage/Desktop/QGRS Mapper_Analysis/P.falciparum_analysis.csv") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if "P" in row["Chromosome"] and not "RC" in row["Chromosome"]:  # if the chromosome is a contig and not RC
                if row["Telomere?"] == "Yes":  # filter out telomeres
                    continue
                else:
                    PQS_dict[str(row["Location"])].append([str(row["Chromosome"]), int(row["PQS Length"])])
            elif "P" in row["Chromosome"] and "RC" in row["Chromosome"]:  # if the chromosome is a contig and RC
                if row["Telomere?"] == "Yes":  # filter out telomeres
                    continue
                else:
                    PQS_dict_RC[str(row["RC_Location Calculated"])].append([str(row["Chromosome"]), int(row["PQS Length"])])
            elif not "RC" in row["Chromosome"]:  # if not a contig and not RC
                if row["Telomere?"] == "Yes":
                    continue
                else:
                    PQS_dict[str(row["Location"])].append([str(row["Chromosome"]), int(row["PQS Length"])])
            elif "RC" in row["Chromosome"]:  # if not a contig and is RC
                if row["Telomere?"] == "Yes":  # filter out telomeres
                    continue
                else:
                    PQS_dict_RC[str(row["RC_Location Calculated"])].append([str(row["Chromosome"]), int(row["PQS Length"])])

cleanPQSfile()

var_gene_lengths = {"lengths": [], "mean": 0, "median": 0, "std_dev": 0}

# get mean and median length of var genes in the species
def calc_var_length():
    """Retrieve the mean, median, and standard deviation of var/pir gene lengths in the genome."""
    for gene in seqs_info:  # loop over all transcripts in the genome
        #if "pir" in gene.description[2].lower() or "rifin" in gene.description[2].lower() or "stevor" in \
                #gene.description[2].lower():  # use this for pir genes
        if "erythrocyte membrane protein 1" in gene.description[2].lower():
            # gene length = end location - start location
            var_gene_length = abs(gene.description[3][1][1] - gene.description[3][1][0])
            var_gene_lengths["lengths"].append(var_gene_length)
    var_gene_lengths.update({"mean": statistics.mean(var_gene_lengths["lengths"])})
    var_gene_lengths.update({"median": statistics.median(var_gene_lengths["lengths"])})
    var_gene_lengths.update({"std_dev": statistics.stdev(var_gene_lengths["lengths"])})

calc_var_length()

random_var_sizes = []  # stores 1,000,000 simulated var/pir gene sizes
def sample_from_normal_distribution():
    """Generate and sample from a distribution of var/pir gene lengths that is normally distributed around the mean,
    and store these lengths in random_var_sizes. Note that a size of 1 million simulations can cause the program to
    take ~3-5 minutes to execute. A size of 100,000 simulations does not drastically change the results, and takes
    significantly less time to execute.
    """
    random_var_length = numpy.random.normal(loc=var_gene_lengths["mean"], scale=var_gene_lengths["std_dev"], size=1000000)
    for length in random_var_length:
        random_var_sizes.append(length)

sample_from_normal_distribution()

''' retrieves the chromosome lengths and chromosome IDs for the species and adds them to lists
    if the chromosome is a contig, the chromosome ID remains the same (eg: "PADLG01_00_01")
    if the chromosome is not a contig, the chromosome ID is the numeric value of the chromosome (eg: "1", "2", "10", etc) '''

chromosome_sizes = []  # stores the size of each chromosome
chromosome_IDs = []  # stores each chromosome ID (eg: "1", "2", "10", "PADLG01_00_01")
def get_chromosome_sizes():
    """Retrieve chromosome lengths and chromosome IDs for the species, and add them to chromosome_sizes and
    chromosome_IDs. If the chromosome number is less than 10, the chromosomeID is the chromosome number without the
    leading "0" (eg: "1"). If the chromosome number is greater than or equal to 10, the chromosomeID is the same as the
    chromosome number (eg: "10"). If the chromosome is a contig, the chromosome ID is the full fasta name (eg: PADLG01_00_01).
    """
    for seq_record in SeqIO.parse(
            "/Users/HunterGage/Desktop/Whole Genome Sequences/PlasmoDB-44_Pfalciparum3d7_Genome.txt", "fasta"):
        length = seq_record.description[seq_record.description.find("length="):]
        length = length[length.find("=") + 1:]
        length = length[:length.find(" ")]
        chromosome = seq_record.id
        chromosome = chromosome[chromosome.find("_") + 1:]
        if "API" in chromosome:  # ignore apicoplast DNA
            continue
        if "M" in chromosome:  # ignore mitochondrial DNA
            continue
        if "00_" in chromosome:  # if the "chromosome" is a contig, chromosomeID should be full fasta name
                                 # (eg: PADLG01_00_01)
            chromosome = seq_record.id
        elif int(chromosome[0]) == 0 and "_" in chromosome:  # if the chromosome number is < 10, chromosomeID should be
                                                             # the chromosome number without leading 0 (eg: "1")
            chromosome = chromosome[1:chromosome.find("_")]
        elif int(chromosome[0]) == 0 and not "_" in chromosome:
            chromosome = chromosome[1]
        elif int(chromosome[0]) == 1 and "_" in chromosome:  # if the chromosome number is >= 10, chromosomeID should be
                                                             # the chromosome number (eg: "10")
            chromosome = chromosome[0:chromosome.find("_")]
        elif int(chromosome[0]) == 1 and not "_" in chromosome:
            chromosome = chromosome[0:]
        chromosome_IDs.append(chromosome)
        chromosome_sizes.append(length)

get_chromosome_sizes()

def calculate_total_genome_size():
    """Add the lengths of all of the chromosomes/contigs to calculate the total size of the genome in bp.

    :return: the total size of the genome in bp
    """
    total_size = 0
    for chromosome_size in chromosome_sizes:
        total_size += int(chromosome_size)
    return total_size

total_size = calculate_total_genome_size()

chromosomes_picked = {}  # track chromosomes picked for testing purposes
def pick_chromosome(var_gene_position):
    """Choose a chromosome based on chromosome lengths, so that larger chromosomes are more likely to be chosen

    :param var_gene_position: a location in the genome, between 0 and the total genome size
    :return: a chromosome/contig number, between 1 and the total number of chromosomes and contigs
    """
    sum = 0
    chromosome = 1
    for chr_size in chromosome_sizes:
        sum += int(chr_size)
        if sum >= var_gene_position:
            if str(chromosome) not in chromosomes_picked:
                chromosomes_picked.update({str(chromosome): 1})
            else:
                chromosomes_picked[str(chromosome)] += 1
            return chromosome
        chromosome += 1

def g4s_on_chromosome(chromosome):
    """Retrieve a list of all of the PQSs on the chosen chromosome. Store PQS location, chromosome, and length in
    dictionaries.
    Dictionary structure: key = the PQS start location
                          value = a list, where index 0 = chromosome ID, and index 1 = PQS length

    :param chromosome: a chromosomeID to find the PQSs on
    :return: a dictionary storing PQS start locations, chromosome IDs, and PQS lengths for all PQSs on the chromosome
    """
    list_of_G4s_on_chromosome = defaultdict(list)
    for g4 in PQS_dict:
        for item in PQS_dict[g4]:
            if str(item[0]) == str(chromosome):
                list_of_G4s_on_chromosome[g4].append([chromosome, item[1]])
    for g4 in PQS_dict_RC:
        for item in PQS_dict_RC[g4]:
            if str(item[0][:item[0].find("_RC")]) == str(chromosome):
                list_of_G4s_on_chromosome[g4].append([chromosome, item[1]])
    return list_of_G4s_on_chromosome

# finds the distance to nearest PQS when given a chromosome, position, and length of a var gene
def find_nearest_g4(chromosome, position, length, start_end):
    """Compute the distance to the nearest PQS when given the chromosome, position, and length of a var/pir gene. If a
    PQS is located within the coding sequence of a pir/var gene, the distance to the nearest PQS = 0. "

    :param chromosome: a chromosomeID that the var/pir gene is located on
    :param position: the location of the var/pir gene on the chromosome
    :param length: the size of the var/pir gene, from the beginning to the end of the coding sequence
    :param start_end: "start" if the "position" refers to the start of the var/pir gene, and "end" if the "position"
                refers to the end of the var/pir gene
    :return: a list, where index 0 = the distance of the var/pir gene from the nearest PQS, and index 1 = the start site
                of the PQS that the var/pir gene is nearest to. The function returns "None" if there are no PQSs on the
                chosen chromosome.
    """
    closest = 1000000000  # track the distance from the nearest PQS
    closest_site = ["no g4 on this chromosome"]  # track the location of the nearest PQS
    g4_data = g4s_on_chromosome(chromosome)  # retrieve the list of PQSs on the chromosome
    g4_start_sites = []  # store PQS start sites
    g4_end_sites = []  # store PQS end sites

    for g4_start_site in g4_data:  # populate the list of PQS start sites for the chosen chromosome
        for item in g4_data[g4_start_site]:
            g4_start_sites.append(g4_start_site)

    for g4_start_site in g4_data:  # populate the list of PQS end sites for the chosen chromosome
        for item in g4_data[g4_start_site]:
            g4_end_sites.append(int(g4_start_site) + int(item[1]))  # PQS end site = PQS start site + PQS length

    if len(g4_start_sites) < 1 or len(g4_end_sites) < 1:  # if there are no PQSs on the chosen chromosome, return None
        return None

    for i in range(len(g4_start_sites)):  # loop over all PQSs on the chromosome to find the closest one
        # if PQS is in var/pir gene, the distance from nearest PQS is 0
        if start_end == "start" and (position < int(g4_start_sites[i]) < (position + abs(length))):
            closest = 0
            closest_site = g4_start_sites[i]
            break

        # check the start site of the PQS to see if it is closest
        if (abs(int(position) - int(g4_start_sites[i]))) < closest:
            closest = abs(int(position) - int(g4_start_sites[i]))
            closest_site = g4_start_sites[i]

        # check the end site of the PQS to see if it is closest
        if (abs(int(position) - int(g4_end_sites[i])) < closest):
            closest = abs(int(position) - int(g4_end_sites[i]))
            closest_site = g4_end_sites[i]

    return [closest, closest_site]

simulated_distances = []  # store distances between simulated var/pir genes and their nearest PQSs

def do_simulation():
    """This is the main function to calculate distances of randomly generated var/pir genes from their nearest PQS. The
    distances are stored in simulated_distances.
    """
    for var_gene_length in random_var_sizes:  # repeat for each simulated var gene
        chromosome_number = pick_chromosome(random.randint(0, total_size))  # pick a chromosome based on genome length
        chromosome_ID = chromosome_IDs[(int(chromosome_number)-1)]  # find the associated chromosome ID
                                                                    # need to subtract 1 because list index starts at 0
        position = random.randint(0, int(chromosome_sizes[int(chromosome_number)-1]))  # choose a random position on the
                                                                                       # chosen chromosome
        if var_gene_length < 0:  # very rarely, the simulation might generate a negative var/pir gene size
                                 # if this occurs, take the "start" position to be the position + length
            nearest_from_start_of_var = find_nearest_g4(chromosome_ID, (position + var_gene_length), abs(var_gene_length), "start")
            if nearest_from_start_of_var is None:  # if no PQSs were found on chromosome, continue to next var/pir gene
                continue
            else:
                nearest_from_start_of_var = nearest_from_start_of_var[0]
            nearest_from_end_of_var = find_nearest_g4(chromosome_ID, position, abs(var_gene_length), "end")[0]
            nearest = min(nearest_from_start_of_var, nearest_from_end_of_var)
            simulated_distances.append(nearest)

        elif var_gene_length >= 0:
            nearest_from_start_of_var = find_nearest_g4(chromosome_ID, position, abs(var_gene_length), "start")
            if nearest_from_start_of_var is None:  # if no PQSs were found on chromosome, continue to next var/pir gene
                continue
            else:
                nearest_from_start_of_var = nearest_from_start_of_var[0]
            nearest_from_end_of_var = find_nearest_g4(chromosome_ID, (position + abs(var_gene_length)), abs(var_gene_length), "end")[0]
            nearest = min(nearest_from_start_of_var, nearest_from_end_of_var)
            #if nearest > 1000000: use if excluding outliers
                #continue
            simulated_distances.append(nearest)

do_simulation()

actual_distances = []  # store distances between actual var/pir genes and their nearest PQSs

def do_my_data():
    """This is the main function to calculate distances of actual var/pir genes from their nearest PQS. The distances
    are stored in actual_distances.

    var_genes dictionary structure: key = the var gene start location
                                    value = a list of dict, with keys/values for "chromosome" and "length"
    """
    var_genes = defaultdict(list)  # store data for var/pir gene location, chromosome and length

    # populate var_genes list with data for var genes
    for seq in seqs_info:
        #if "pir" in seq.description[2].lower() or "rifin" in seq.description[2].lower() or "stevor" in \
                #seq.description[2].lower(): # use for pir genes
        if "erythrocyte membrane protein 1" in seq.description[2].lower():
            chromosome = seq.description[3][0][seq.description[3][0].find("=")+1:]  # get chromosome var/pir gene is on

            if "_00_" in chromosome:  # if chromosome is a contig, leave chromosomeID as is (eg: "PADLG01_00_01")
                chromosome = chromosome
                var_genes[seq.description[3][1][0]].append([{"chromosome": chromosome,
                                                             "length": abs(
                                                                 seq.description[3][1][1] - seq.description[3][1][0])}])
                continue
            chromosome = chromosome[chromosome.find("_")+1:]  # if chromosome is not a contig, represent the ID as
                                                              # a number (eg: "1" or "11")
            if int(chromosome[0]) == 0 and "_" in chromosome:  # if chromosome number is less than 10, drop leading "0"
                chromosome = chromosome[1:chromosome.find("_")]
            elif int(chromosome[0]) == 0 and not "_" in chromosome:
                chromosome = chromosome[1]
            elif int(chromosome[0]) == 1 and "_" in chromosome:  # if chromosome number is >= 10, leave as is
                chromosome = chromosome[0:chromosome.find("_")]
            elif int(chromosome[0]) == 1 and not "_" in chromosome:
                chromosome = chromosome[0:]

            var_genes[seq.description[3][1][0]].append([{"chromosome": chromosome,
                                           "length": abs(seq.description[3][1][1] - seq.description[3][1][0])}])

    # find distance from each var/pir gene to the nearest PQS
    for var_gene in var_genes:
        for item in var_genes[var_gene]:
            nearest_from_start_of_var = find_nearest_g4(item[0]["chromosome"], var_gene, item[0]["length"], "start")
            if nearest_from_start_of_var is None:  # if no PQSs on the chromosome, skip the gene
                print("no PQSs found on chr" + item[0]["chromosome"])
                continue
            else:
                nearest_from_start_of_var = nearest_from_start_of_var[0]
            nearest_from_end_of_var = find_nearest_g4(item[0]["chromosome"], (int(var_gene) + int(item[0]["length"])),
                                                      item[0]["length"], "end")[0]
            nearest = min(nearest_from_start_of_var, nearest_from_end_of_var)
            #if nearest > 1000000:  use if excluding outliers
                #continue
            actual_distances.append(nearest)

do_my_data()

def get_confidence_interval(dataset):
    """Calculates a 95% confidence interval around the median of a list.

    :param dataset: a list of distances between var/pir genes and their nearest PQSs
    :return: a list, where index 0 = the lower bound of the 95% CI, and index 1 = the upper bound of the 95% CI
    """
    ordered_dataset = sorted(dataset)
    sample_size = len(ordered_dataset)
    lower_95_limit_index = round((sample_size / 2) - ((1.96 * (sample_size **0.5)) / 2)) - 1
    upper_95_limit_index = round(1 + (sample_size / 2) + ((1.96 * (sample_size **0.5)) / 2)) - 1
    lower_95_limit_value = ordered_dataset[lower_95_limit_index]
    upper_95_limit_value = ordered_dataset[upper_95_limit_index]
    return [lower_95_limit_value, upper_95_limit_value]

def display_distributions():
    # plot the distribution of simulated var/pir gene lengths, normally distributed around the mean
    plt.hist(random_var_sizes, bins=300)
    plt.xlabel("var/pir gene lengths")
    plt.ylabel("count")
    plt.show()

    chromosomes_picked_to_sort = []
    for chr in chromosomes_picked:
        chromosomes_picked_to_sort.append(int(chr))
    chromosomes_picked_to_sort = sorted(chromosomes_picked_to_sort)
    summary_dict = {}
    for chromosome in chromosomes_picked_to_sort:
        summary_dict.update({chromosome: 0})
    for chr in chromosomes_picked:
        summary_dict[int(chr)] = chromosomes_picked[chr]

    # plot the distribution of chromosomes chosen; larger chromosomes should be chosen more frequently
    plt.bar(list(summary_dict.keys()), list(summary_dict.values()))
    plt.ylabel("count")
    plt.xlabel("chromosome")
    plt.show()

    # plot the distribution of distances between var/pir genes and their nearest PQS, for the simulated dataset
    plt.hist(simulated_distances, bins=150)
    plt.xlabel("simulated distance from nearest PQS (bp)")
    plt.ylabel("count")
    plt.show()

    # plot the distribution of distances between var/pir genes and their nearest PQS, for the actual dataset
    plt.hist(actual_distances, bins=150)
    plt.xlabel("actual distance from nearest PQS (bp)")
    plt.ylabel("count")
    plt.show()


# print the results to the console
print("mean distance in null dataset: ", statistics.mean(simulated_distances))
print("mean distance in actual dataset: ", statistics.mean(actual_distances))
print("median distance in null dataset: ", statistics.median(simulated_distances))
print("median distance in actual dataset: ", statistics.median(actual_distances))
print("95% confidence interval for actual dataset: ", get_confidence_interval(actual_distances)[0], ",",
      get_confidence_interval(actual_distances)[1])
print(stats.ttest_ind(simulated_distances, actual_distances, equal_var=False))
print("no. genes in family: ", len(var_gene_lengths["lengths"]))
print("no. genes in family on a chromosome with at least one PQS: ", len(actual_distances))
display_distributions()
