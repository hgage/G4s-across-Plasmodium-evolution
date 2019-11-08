"""This script takes the data from FindAnnotation_QGRSMapper.py and updates the "ups" / "down" and "sense" / "antisense"
 designations. The PQS and nearest gene are considered "sense" if they are on the same strand, and "antisense" if
 they are on opposite strands. The data is written to a new .csv file."""

import csv

with open("/Users/HunterGage/PycharmProjects/9/Practice/P.chabaudi_QGRSannotation.csv", 'r') as csvfile:
    reader = csv.reader(csvfile)
    rows = []
    for row in reader:
        if row[0] == "Location":
            continue
        if row[2] == "reverse-complement":  # PQS is on complement strand
            if row[5] == "+":  # gene is on sense strand
                row[5] = "Antisense"  # PQS and gene are on different strands
                if row[3].lower() == "ups":
                    row[3] = "Ups"
                    if not int(row[6]) > 0:
                        row[6] = (-1) * int(row[6])
                elif row[3].lower() == "down":
                    row[3] = "Down"
                    if not int(row[6]) < 0:
                        row[6] = (-1) * int(row[6])
            elif row[5] == "-":  # gene is on antisense strand
                row[5] = "Sense"  # PQS and gene are both on the antisense strand
                if row[3].lower() == "ups":
                    row[3] = "Down"
                    if not int(row[6]) < 0:
                        row[6] = (-1) * int(row[6])
                elif row[3].lower() == "down":
                    row[3] = "Ups"
                    if not int(row[6]) > 0:
                        row[6] = (-1) * int(row[6])
        elif row[2] == "sense":  # PQS is on sense strand
            if row[5] == "+":  # gene is on sense strand
                row[5] = "Sense"  # PQS and gene are both on sense strand
                if row[3].lower() == "ups":
                    row[3] = "Ups"
                    if not int(row[6]) > 0:
                        row[6] = (-1) * int(row[6])
                elif row[3].lower() == "down":
                    row[3] = "Down"
                    if not int(row[6]) < 0:
                        row[6] = (-1) * int(row[6])
            elif row[5] == "-":  # gene is on antisense strand
                row[5] = "Antisense"  # PQS and gene are on different strands
                if row[3].lower() == "ups":
                    row[3] = "Down"
                    if not int(row[6]) < 0:
                        row[6] = (-1) * int(row[6])
                elif row[3].lower() == "down":
                    row[3] = "Ups"
                    if not int(row[6]) > 0:
                        row[6] = (-1) * int(row[6])
        rows.append(row)
    print(rows)

# write the data to a new .csv file
with open('P.chabaudi_QGRScleanup.csv', 'w') as summaryfile:
    headers = ["Location", "Chromosome", "Sense or reverse-complement", "Coding, upstream, or downstream?",
               "Nearest gene ID", "Strand", "Distance", "Gene Annotation"]
    filewriter = csv.writer(summaryfile)
    filewriter.writerow(headers)
    for data in rows:
        filewriter.writerow(data)
