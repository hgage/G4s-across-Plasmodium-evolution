"""This script takes a data set of PQS locations and parses through the annotated coding sequences to find the nearest
gene to each PQS. It calculates the distance of each PQS to its nearest gene and writes the data to a csv file."""

# import some useful libraries
from Bio import SeqIO
import csv
from collections import defaultdict

# list of sequence record objects for each gene annotation in the fasta file
seqs_info = []
for seq_record in SeqIO.parse("/Users/HunterGage/Desktop/Gene Annotations/PlasmoDB-44_Pfalciparum3d7_AnnotatedCDSs.txt", "fasta"):
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
PQS_dict_RC = defaultdict(list)  # stores data for PQSs on the reverse complement strand

def cleanPQSfile():
    """Read from csv file containing PQSs and their locations. Store PQS location, chromosome, and length in PQS_dict
        and PQS_dict_RC.
        Dictionary structure: key = PQS start location
                              value = a list, where index 0 = chromosome number, and index 1 = PQS length
    """
    with open("/Users/HunterGage/Desktop/QGRS Mapper_Analysis/P.falciparum_analysis.csv") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if "P" in row["Chromosome"] and not "RC" in row["Chromosome"]: # if the chromosome is a contig and not RC
                PQS_dict[str(row["Location"])].append([str(row["Chromosome"])])
            elif "P" in row["Chromosome"] and "RC" in row["Chromosome"]:  # if the chromosome is a contig and RC
                PQS_dict_RC[str(row["RC_Location Calculated"])].append([str(row["Chromosome"])])
            elif not "RC" in row["Chromosome"]:  # if the chromosome is not a contig and not RC
                PQS_dict[str(row["Location"])].append([int(row["Chromosome"])])
            elif "RC" in row["Chromosome"]:  # if the chromosome is not a contig and is RC
                PQS_dict_RC[str(row["RC_Location Calculated"])].append([str(row["Chromosome"])])

cleanPQSfile()

# adds identifier to each chromosome in the PQS dict so that they can be compared to the chromosomes in the FASTA file
clean_PQS_dict = defaultdict(list)
clean_PQS_dict_RC = defaultdict(list)
def cleanchromosomes():
    """Add identifier to each chromosome in the PQS dicts so that they can be compared to the chromosomes in FASTA file.
    For example, chromosome "1" must become "PADLG01_01" to follow the naming conventions in the FASTA file.
    """
    for location in PQS_dict:  # loop over all PQSs on the sense strand
        for entry in PQS_dict[location]:
            if "P" in str(entry[0]) and not "RC" in str(entry[0]):  # if the chromosome is a contig and not RC
                clean_PQS_dict[location].append([entry[0]])
            elif entry[0] < 10:
                clean_PQS_dict[location].append(["Pf3D7_0" + str(entry[0]) + "_v3"])  # add identifier to match format
                                                                                      # of chromosome header in PlasmoDB
            elif entry[0] >= 10:
                clean_PQS_dict[location].append(["Pf3D7_" + str(entry[0]) + "_v3"])
    for location in PQS_dict_RC:  # loop over all PQSs on the reverse complement strand
        for entry in PQS_dict_RC[location]:
            if "P" in str(entry[0]) and "RC" in str(entry[0]):
                clean_PQS_dict_RC[location].append([entry[0][:entry[0].find("_RC")]])
            elif int(entry[0][:entry[0].find("_RC")]) < 10:
                clean_PQS_dict_RC[location].append(["Pf3D7_0" + str(entry[0][:entry[0].find("_RC")]) + "_v3"])
            elif int(entry[0][:entry[0].find("_RC")]) >= 10:
                clean_PQS_dict_RC[location].append(["Pf3D7_" + str(entry[0][:entry[0].find("_RC")]) + "_v3"])

cleanchromosomes()

summary_data = defaultdict(list)  # store data for nearest gene to each PQS on the sense strand
def parseannotations():
    """Parse the FASTA file for each PQS and retrieve information for chromosome, position relative to nearest gene,
    nearest gene ID, strand, distance from nearest gene, and nearest gene annotation. Add the data to summary_data.
    Dictionary structure:     key = PQS start location
                              value = a list, where index 0 = a dict, with keys: "Chromosome", "Coding, upstream,
                                   or downstream", "Nearest gene ID", "Strand", "Distance", and "Gene annotation"
    """

    #looping over all of the PQS that we wish to search
    for PQS_location in clean_PQS_dict:
        for entry in clean_PQS_dict[PQS_location]:
            index = 0
            seqs_in_chromosome = []
            for seq in seqs_info:
                if entry[0] == seq.description[3][0][seq.description[3][0].find("=")+1:]:
                    seqs_in_chromosome.append(seq)

            #if the chromosome / contig doesn't have any annotated genes
            if len(seqs_in_chromosome) == 0:
                print("warning: " + entry[0] + " doesn't have any annotated genes")
                summary_data[str(PQS_location)].append([{
                    "Chromosome": entry[0],
                    "Coding, upstream, or downstream?": "None",
                    "Nearest gene ID": "None",
                    "Strand": "None",
                    "Distance": "None",
                    "Gene annotation": "None"}])
                continue

            #looping over all of the annotations in the FASTA file from plasmodb that are on the same chromosome as the PQS
            for seq in seqs_in_chromosome:
                try:
                    # if PQS comes before any annotated genes in the chromosome
                    if int(PQS_location) <= seq.description[3][1][0]:
                        summary_data[str(PQS_location)].append([{
                            "Chromosome": seq.description[3][0][seq.description[3][0].find("=") + 1:],
                            "Coding, upstream, or downstream?": "Ups",
                            "Nearest gene ID": seq.description[0],
                            "Strand": seq.description[3][2][:seq.description[3][2].find("\)")],
                            "Distance": seq.description[3][1][0] - int(PQS_location),
                            "Gene annotation": seq.description[2][seq.description[2].find("=") + 1:]}])
                        break

                    #if the location of the PQS is in a gene
                    elif seq.description[3][1][0] <= int(PQS_location) <= seq.description[3][1][1]:
                        summary_data[str(PQS_location)].append([{"Chromosome": seq.description[3][0][seq.description[3][0].find("=")+1:],
                                                                 "Coding, upstream, or downstream?": "In",
                                                                 "Nearest gene ID": seq.description[0],
                                                                 "Strand": seq.description[3][2][:seq.description[3][2].find("\)")],
                                                                 "Distance": 0,
                                                                 "Gene annotation": seq.description[2][
                                                                                    seq.description[2].find("=") + 1:]}])
                        break

                    #if the location of the PQS is between genes
                    elif int(seq.description[3][1][1]) <= int(PQS_location) <= int(seqs_in_chromosome[index+1].description[3][1][0]):
                        distance_from_first_gene = int(PQS_location) - int(seq.description[3][1][1])
                        distance_from_second_gene = int(seqs_in_chromosome[index+1].description[3][1][0]) - int(PQS_location)
                        # if PQS is closer to the first gene than the second gene
                        if distance_from_first_gene <= distance_from_second_gene:
                            summary_data[str(PQS_location)].append([{
                                "Chromosome": seq.description[3][0][seq.description[3][0].find("=") + 1:],
                                "Coding, upstream, or downstream?": "Down",
                                "Nearest gene ID": seq.description[0],
                                "Strand": seq.description[3][2][:seq.description[3][2].find("\)")],
                                "Distance": (-1) * distance_from_first_gene,
                                "Gene annotation": seq.description[2][seq.description[2].find("=") + 1:]}])
                            break
                        # if PQS is closer to the second gene than the first gene
                        elif distance_from_first_gene > distance_from_second_gene:
                            summary_data[str(PQS_location)].append([{
                                "Chromosome": seqs_in_chromosome[index+1].description[3][0][seq.description[3][0].find("=") + 1:],
                                "Coding, upstream, or downstream?": "Ups",
                                "Nearest gene ID": seqs_in_chromosome[index+1].description[0],
                                "Strand": seqs_in_chromosome[index+1].description[3][2][:seq.description[3][2].find("\)")],
                                "Distance": distance_from_second_gene,
                                "Gene annotation": seqs_in_chromosome[index + 1].description[2][
                                                   seq.description[2].find("=") + 1:]}])
                            break

                    index += 1

                except:
                    # if PQS is after all annotated genes in the chromosome
                    summary_data[str(PQS_location)].append([{
                        "Chromosome": seq.description[3][0][seq.description[3][0].find("=") + 1:],
                        "Coding, upstream, or downstream?": "Down",
                        "Nearest gene ID": seq.description[0],
                        "Strand": seq.description[3][2][:seq.description[3][2].find("\)")],
                        "Distance": seq.description[3][1][1] - int(PQS_location),
                        "Gene annotation": seq.description[2][seq.description[2].find("=") + 1:]}])
                    break

    return summary_data

parseannotations()

summary_data_RC = defaultdict(list)  # store data for nearest gene to each PQS on the reverse complement strand
def parseannotations_RC():
    """Parse the FASTA file for each PQS and retrieve information for chromosome, position relative to nearest gene,
        nearest gene ID, strand, distance from nearest gene, and nearest gene annotation. Add the data to summary_data_RC.
        Dictionary structure:     key = PQS start location
                                  value = a list, where index 0 = a dict, with keys: "Chromosome", "Coding, upstream,
                                       or downstream", "Nearest gene ID", "Strand", "Distance", and "Gene annotation"
    """
    #looping over all of the PQS that we wish to search
    for PQS_location in clean_PQS_dict_RC:
        for entry in clean_PQS_dict_RC[PQS_location]:
            index = 0
            seqs_in_chromosome = []
            for seq in seqs_info:
                if entry[0] == seq.description[3][0][seq.description[3][0].find("=")+1:]:
                    seqs_in_chromosome.append(seq)

            # if the chromosome / contig doesn't have any annotated genes
            if len(seqs_in_chromosome) == 0:
                print("warning: " + entry[0] + " doesn't have any annotated genes")
                summary_data_RC[str(PQS_location)].append([{
                    "Chromosome": entry[0],
                    "Coding, upstream, or downstream?": "None",
                    "Nearest gene ID": "None",
                    "Strand": "None",
                    "Distance": "None",
                    "Gene annotation": "None"}])
                continue

            # looping over all of the annotations in the FASTA file from plasmodb that are on the same chromosome as the PQS
            for seq in seqs_in_chromosome:
                try:
                    # if PQS comes before any annotated genes in the chromosome
                    if int(PQS_location) <= seq.description[3][1][0]:
                        summary_data_RC[str(PQS_location)].append([{
                            "Chromosome": seq.description[3][0][seq.description[3][0].find("=") + 1:],
                            "Coding, upstream, or downstream?": "Ups",
                            "Nearest gene ID": seq.description[0],
                            "Strand": seq.description[3][2][:seq.description[3][2].find("\)")],
                            "Distance": seq.description[3][1][0] - int(PQS_location),
                            "Gene annotation": seq.description[2][seq.description[2].find("=") + 1:]}])
                        break

                    #if the location of the PQS is in a gene
                    elif seq.description[3][1][0] <= int(PQS_location) <= seq.description[3][1][1]:
                        summary_data_RC[str(PQS_location)].append([{"Chromosome": seq.description[3][0][seq.description[3][0].find("=")+1:],
                                                                 "Coding, upstream, or downstream?": "In",
                                                                 "Nearest gene ID": seq.description[0],
                                                                 "Strand": seq.description[3][2][:seq.description[3][2].find("\)")],
                                                                 "Distance": 0,
                                                                 "Gene annotation": seq.description[2][
                                                                                    seq.description[2].find("=") + 1:]}])
                        break

                    #if the location of the PQS is between genes
                    elif int(seq.description[3][1][1]) <= int(PQS_location) <= int(seqs_in_chromosome[index+1].description[3][1][0]):
                        distance_from_first_gene = int(PQS_location) - int(seq.description[3][1][1])
                        distance_from_second_gene = int(seqs_in_chromosome[index+1].description[3][1][0]) - int(PQS_location)
                        # if PQS is closer to the first gene than the second gene
                        if distance_from_first_gene <= distance_from_second_gene:
                            summary_data_RC[str(PQS_location)].append([{
                                "Chromosome": seq.description[3][0][seq.description[3][0].find("=") + 1:],
                                "Coding, upstream, or downstream?": "Down",
                                "Nearest gene ID": seq.description[0],
                                "Strand": seq.description[3][2][:seq.description[3][2].find("\)")],
                                "Distance": (-1) * distance_from_first_gene,
                                "Gene annotation": seq.description[2][seq.description[2].find("=") + 1:]}])
                            break
                        # if PQS is closer to the second gene than the first gene
                        elif distance_from_first_gene > distance_from_second_gene:
                            summary_data_RC[str(PQS_location)].append([{
                                "Chromosome": seqs_in_chromosome[index+1].description[3][0][seq.description[3][0].find("=") + 1:],
                                "Coding, upstream, or downstream?": "Ups",
                                "Nearest gene ID": seqs_in_chromosome[index+1].description[0],
                                "Strand": seqs_in_chromosome[index+1].description[3][2][:seq.description[3][2].find("\)")],
                                "Distance": distance_from_second_gene,
                                "Gene annotation": seqs_in_chromosome[index + 1].description[2][
                                                   seq.description[2].find("=") + 1:]}])
                            break

                    index += 1

                except:
                    # if PQS is after all annotated genes in the chromosome
                    summary_data_RC[str(PQS_location)].append([{
                        "Chromosome": seq.description[3][0][seq.description[3][0].find("=") + 1:],
                        "Coding, upstream, or downstream?": "Down",
                        "Nearest gene ID": seq.description[0],
                        "Strand": seq.description[3][2][:seq.description[3][2].find("\)")],
                        "Distance": seq.description[3][1][1] - int(PQS_location),
                        "Gene annotation": seq.description[2][seq.description[2].find("=") + 1:]}])
                    break

    return summary_data_RC

parseannotations_RC()

#creates csv file of summary data
def makecsv():
    """Create csv file of the data contained in summary_data and summary_data_RC, with columns for: "Location",
    "Chromosome", "Sense or reverse-complement", "Coding, upstream, or downstream?", "Nearest gene ID", "Strand",
    "Distance", and "Gene Annotation"
    """
    d = []  # store the data compiled from summary_data and summary_data_RC, to be written to csv
    for data in summary_data:  # add the data for the PQSs on the sense strand
        for item in summary_data[data]:
            d.append(
                [data, item[0]["Chromosome"], "sense", item[0]["Coding, upstream, or downstream?"],
                 item[0]["Nearest gene ID"], item[0]["Strand"],
                 item[0]["Distance"], item[0]["Gene annotation"]])

    for data in summary_data_RC:  # add the data for the PQSs on the reverse complement strand
        for item in summary_data_RC[data]:
            d.append(
                [data, item[0]["Chromosome"], "reverse-complement", item[0]["Coding, upstream, or downstream?"],
                 item[0]["Nearest gene ID"], item[0]["Strand"],
                 item[0]["Distance"], item[0]["Gene annotation"]])

    # write the data to csv file
    with open ('P.falciparum_QGRSannotation.csv', 'w') as summaryfile:
        headers = ["Location", "Chromosome", "Sense or reverse-complement", "Coding, upstream, or downstream?", "Nearest gene ID", "Strand", "Distance", "Gene Annotation"]
        filewriter = csv.writer(summaryfile)
        filewriter.writerow(headers)
        for data in d:
            filewriter.writerow(data)

makecsv()