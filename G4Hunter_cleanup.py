"""This script takes the data from FindAnnotation_G4Hunter.py and updates the "ups" / "down" and "sense" / "antisense"
designations. The PQS and nearest gene are considered "sense if they are on the same strand, and "antisense" if they
are on opposite strands. The data is written to a new .csv file."""

import csv

with open("/Users/HunterGage/Desktop/G4 Hunter_Analysis/P.chabaudi_G4analysis.csv", 'r') as csvfile:
    reader = csv.reader(csvfile)
    rows = []
    for row in reader:
        if row[0] == "Chromosome":
            continue
        if float(row[3]) < 0:  # PQS is on complement strand
            if row[8] == "+":  # gene is on sense strand
                row[8] = "Antisense"  # PQS and gene are on different strands
                if row[6].lower() == "ups":
                    row[6] = "Ups"
                    if not int(row[9]) > 0:
                        row[9] = (-1) * int(row[9])
                elif row[6].lower() == "down":
                    row[6] = "Down"
                    if not int(row[9]) < 0:
                        row[9] = (-1) * int(row[9])
            elif row[8] == "-":  # gene is on antisense strand
                row[8] = "Sense"  # PQS and gene are both on the antisense strand
                if row[6].lower() == "ups":
                    row[6] = "Down"
                    if not int(row[9]) < 0:
                        row[9] = (-1) * int(row[9])
                elif row[6].lower() == "down":
                    row[6] = "Ups"
                    if not int(row[9]) > 0:
                        row[9] = (-1) * int(row[9])
        elif float(row[3]) > 0:  # PQS is on sense strand
            if row[8] == "+":  # gene is on sense strand
                row[8] = "Sense"  # PQS and gene are both on sense strand
                if row[6].lower() == "ups":
                    row[6] = "Ups"
                    if not int(row[9]) > 0:
                        row[9] = (-1) * int(row[9])
                elif row[6].lower() == "down":
                    row[6] = "Down"
                    if not int(row[9]) < 0:
                        row[9] = (-1) * int(row[9])
            elif row[8] == "-":  # gene is on antisense strand
                row[8] = "Antisense"  # PQS and gene are on different strands
                if row[6].lower() == "ups":
                    row[6] = "Down"
                    if not int(row[9]) < 0:
                        row[9] = (-1) * int(row[9])
                elif row[6].lower() == "down":
                    row[6] = "Ups"
                    if not int(row[9]) > 0:
                        row[9] = (-1) * int(row[9])
        rows.append(row)
    print(rows)

# write the data to a new.csv file
with open('P.chabaudi_G4cleanup.csv', 'w') as summaryfile:
    headers = ["Chromosome", "Location", "Length", "Score", "Abs Score", "Sequence", "Coding, upstream, or downstream?", "Nearest gene ID", "Strand", "Distance",
               "Gene Annotation"]
    filewriter = csv.writer(summaryfile)
    filewriter.writerow(headers)
    for data in rows:
        filewriter.writerow(data)