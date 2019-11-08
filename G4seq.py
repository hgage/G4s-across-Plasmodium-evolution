"""This script compares the concordance of QGRS Mapper and experimental G4-seq data by searching for instances where
 the locations of PQSs predicted by QGRS Mapper overlap with a G4-seq window."""

# import some useful libraries
import csv
from collections import defaultdict

PQS_dict = defaultdict(list)  # stores PQSs on the sense strand
PQS_dict_RC = defaultdict(list)  # stores PQSs on the reverse complement strand
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
                PQS_dict[str(row["Location"])].append([str(row["Chromosome"]), int(row["PQS Length"])])
            elif "P" in row["Chromosome"] and "RC" in row["Chromosome"]:  # if the chromosome is a contig and RC
                PQS_dict_RC[str(row["RC_Location Calculated"])].append([str(row["Chromosome"]), int(row["PQS Length"])])
            elif not "RC" in row["Chromosome"]:  # if not a contig and not RC
                PQS_dict[str(row["Location"])].append([str(row["Chromosome"]), int(row["PQS Length"])])
            elif "RC" in row["Chromosome"]:  # if not a contig and RC
                PQS_dict_RC[str(row["RC_Location Calculated"])].append([str(row["Chromosome"]), int(row["PQS Length"])])

cleanPQSfile()

G4seqs = []  # stores G4 seq windows on the sense strand
G4seqs_reverse = []  # stores G4 seq windows on the reverse complement strand
def cleanG4seqfile():
    """Read from csv files containing G4 seq windows and their locations. Store chromosome, start location, and end
    location in G4seqs and G4seqs_reverse.
    """
    with open("/Users/HunterGage/Desktop/Whole Genome Sequences/G4seq_forward.csv") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            data = {"chromosome": row["Chromosome"], "start": row["Start"], "end": row["End"]}
            G4seqs.append(data)
    with open("/Users/HunterGage/Desktop/Whole Genome Sequences/G4seq_reverse.csv") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            data = {"chromosome": row["Chromosome"], "start": row["Start"], "end": row["End"]}
            G4seqs_reverse.append(data)
cleanG4seqfile()

data = []  # for each PQS, stores whether the PQS was found, as well as the G4-seq window it was found within
def compare_lists():
    """For each PQS found by QGRS Mapper, determine whether its location overlaps with a G4-seq window, and store this
    information in data. If the PQS overlaps with a G4-seq window, also retrieve the location of G4-seq window.
    """
    for PQS in PQS_dict:  # loop over all PQSs on the sense strand
        for item in PQS_dict[PQS]:
            found = False  # track whether the PQS was found within a G4-seq window
            for G4_seq in G4seqs:  # loop over all of the G4-seq windows on the sense strand
                if item[0] == G4_seq["chromosome"]:  # only consider PQSs and G4-seq windows on the same chromosome
                    if int(G4_seq["start"]) <= int(PQS) <= int(G4_seq["end"]):
                        data.append(str(PQS) + ", " + str(item[0]) + " found " + str(G4_seq["start"]) + "-"
                                    + str(G4_seq["end"]))
                        found = True
                        break
                    elif int(G4_seq["start"]) <= (int(PQS) + int(item[1])) <= int(G4_seq["end"]):
                        data.append(str(PQS) + ", " + str(item[0]) + " found " + str(G4_seq["start"]) + "-"
                                    + str(G4_seq["end"]))
                        found = True
                        break
            if found == False:
                data.append(str(PQS) + ", " + str(item[0]) + " not found")
    for PQS in PQS_dict_RC:  # loop over all PQSs on the reverse complement strand
        for item in PQS_dict_RC[PQS]:
            found = False  # track whether the PQS was found within a G4-seq window
            for G4_seq in G4seqs_reverse:  # loop over all of the G4-seq windows on the reverse complement strand
                if item[0] == G4_seq["chromosome"]:  # only consider PQSs and G4-seq windows on the same chromosome
                    if int(G4_seq["start"]) <= int(PQS) <= int(G4_seq["end"]):
                        data.append(str(PQS) + ", " + str(item[0]) + " found " + str(G4_seq["start"]) + "-"
                                    + str(G4_seq["end"]))
                        found = True
                        break
                    elif int(G4_seq["start"]) <= (int(PQS) + int(item[1])) <= int(G4_seq["end"]):
                        data.append(str(PQS) + ", " + str(item[0]) + " found " + str(G4_seq["start"]) + "-"
                                    + str(G4_seq["end"]))
                        found = True
                        break
            if found == False:
                data.append(str(PQS) + ", " + str(item[0]) + " not found")
compare_lists()

G4seq_found_by_algorithm = []  # for each G4-seq window, stores whether it overlaps with a PQS
def compare_lists_G4seq():
    """For each G4-seq window, determine whether its location overlaps with a PQS predicted by QGRS Mapper, and store
    this information in G4seq_found_by_algorithm.
    """
    for G4_seq in G4seqs:  # loop over all of the G4 seq windows on the sense strand
        found = False  # track whether a PQS was found within the G4-seq window
        for PQS in PQS_dict:  # loop over all of the PQSs on the sense strand
            for item in PQS_dict[PQS]:
                if item[0] == G4_seq["chromosome"]:  # only consider G4-seq windows and PQSs on the same chromosome
                    if int(G4_seq["start"]) <= int(PQS) <= int(G4_seq["end"]):
                        G4seq_found_by_algorithm.append(G4_seq)
                        found = True
                        break
                    elif int(G4_seq["start"]) <= (int(PQS) + int(item[1])) <= int(G4_seq["end"]):
                        G4seq_found_by_algorithm.append(G4_seq)
                        found = True
                        break
            if found == True:
                break
    for G4_seq in G4seqs_reverse:  # loop over all of the G4 seq windows on the reverse complement strand
        found = False  # track whether a PQS was found within the G4-seq window
        for PQS in PQS_dict_RC:  # loop over all of the PQSs on the reverse complement strand
            for item in PQS_dict_RC[PQS]:
                if item[0] == G4_seq["chromosome"]:  # only consider G4-seq windows and PQSs on the same chromosome
                    if int(G4_seq["start"]) <= int(PQS) <= int(G4_seq["end"]):
                        G4seq_found_by_algorithm.append(str(PQS) + ", " + str(item[0]) + " found "
                                                        + str(G4_seq["start"]) + "-" + str(G4_seq["end"]))
                        found = True
                        break
                    elif int(G4_seq["start"]) <= (int(PQS) + int(item[1])) <= int(G4_seq["end"]):
                        G4seq_found_by_algorithm.append(str(PQS) + ", " + str(item[0]) + " found "
                                                        + str(G4_seq["start"]) + "-" + str(G4_seq["end"]))
                        found = True
                        break
            if found == True:
                break
compare_lists_G4seq()

for entry in data:  # print the list of PQSs and whether or not they were found within a G4-seq window
    print(entry)

count_G4seqs = 0  # total number of G4-seq windows on the sense strand
for seq in G4seqs:
    count_G4seqs += 1

count_G4seqs_reverse = 0  # total number of G4-seq windows on the reverse complement strand
for seq in G4seqs_reverse:
    count_G4seqs_reverse += 1

print("G4seqs", count_G4seqs, "G4seqs_reverse", count_G4seqs_reverse)

count_G4seq_found_by_algorithm = 0  # total number of G4-seq windows that overlapped with a PQS
for seq in G4seq_found_by_algorithm:
    count_G4seq_found_by_algorithm += 1
print("total G4seq found by algorithm", count_G4seq_found_by_algorithm)


