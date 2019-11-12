"""This script calculates the distance of each PQS from the start codon of its nearest gene, and writes the data to a
csv file."""

# import some useful libraries
import csv
from Bio import SeqIO
from collections import defaultdict

''' store all PQSs in a dictionary, where key = PQS location, and value = a list, where index 0 = nearest gene ID, and
index 1 = nearest gene annotation.'''
G4_locations = defaultdict(list)

with open("/Users/HunterGage/Desktop/G4 Hunter_Analysis/P.praefalciparum_G4analysis.csv", 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if row["Coding, upstream, or downstream?"].lower() == "none":
            G4_locations[str(row["Location"])].append(["None", "None"])
        else:
            G4_locations[str(row["Location"])].append([row["Nearest gene ID"], row["Gene Annotation"]])

seqs_info = []  # list of sequence record objects for each gene annotation in the fasta file
for seq_record in SeqIO.parse("/Users/HunterGage/Desktop/Gene Annotations/PlasmoDB-44_PpraefalciparumG01_AnnotatedCDSs.txt", "fasta"):
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
    """Split the location in the description into a list so that the numeric value can be accessed
            to access start location: seq.description[3][1][0]
            to access end location: seq.description[3][1][1]
            to access strand: seq.description[3][2][:seq.description[3][2].find("\)")]
    """
    for seq in seqs_info:
        seq.description[3] = seq.description[3].replace("(", ":")
        seq.description[3] = seq.description[3].split(":")
        seq.description[3][1] = seq.description[3][1].split("-")
        seq.description[3][1] = list(map(int, seq.description[3][1]))
    return seqs_info

cleanlocations()

summary_data = defaultdict(list)  # store data for the position of each PQS relative to the start codon of nearest gene
def find_position_in_gene():
    """Determine the position of each PQs relative to the start codon of its nearest gene, and store the data in
    summary_data.
    """
    for G4 in G4_locations:  # loop over all PQSs
        for entry in G4_locations[G4]:
            if entry[0] == "None":
                summary_data[G4].append([{"Nearest gene ID": entry[0],
                                          "Position in gene": "None",
                                          "Gene annotation": "None"}])
            else:
                for seq in seqs_info:  # loop over all coding sequences
                    if entry[0] == seq.description[0]:  # nearest gene ID of PQS matches the gene ID of coding sequence
                        if seq.description[3][2][:seq.description[3][2].find("\)")] == "+":  # gene on sense strand
                            position_in_gene = int(G4) - int(seq.description[3][1][0])
                            summary_data[G4].append([{"Nearest gene ID": entry[0],
                                                      "Position in gene": position_in_gene,
                                                      "Gene annotation": entry[1]}])
                        elif seq.description[3][2][:seq.description[3][2].find("\)")] == "-":  # gene on RC strand
                            position_in_gene = int(seq.description[3][1][1]) - int(G4)
                            summary_data[G4].append([{"Nearest gene ID": entry[0],
                                                      "Position in gene": position_in_gene,
                                                      "Gene annotation": entry[1]}])
find_position_in_gene()

def makecsv():
    """Create a csv file from summary_data, with columns for "Location", "Nearest gene ID", "Position in gene", and
    "Gene Annotation".
    """
    d = []  # store data from summary_data, to write to csv
    for data in summary_data:
        for entry in summary_data[data]:
            d.append(
                [data, entry[0]["Nearest gene ID"], entry[0]["Position in gene"], entry[0]["Gene annotation"]])

    # write to csv file
    with open ("/Users/HunterGage/Desktop/G4 Positions Relative to Start Codon/P.praefalciparum_G4Hunter_codon.csv", 'w') as summaryfile:
        headers = ["Location", "Nearest gene ID", "Position in gene", "Gene Annotation"]
        filewriter = csv.writer(summaryfile)
        filewriter.writerow(headers)
        for data in d:
            filewriter.writerow(data)

makecsv()
