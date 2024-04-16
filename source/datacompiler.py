filename = "gene_result.txt"

file = open(filename, "r", encoding="utf-8")
output_file = open("gene_output.csv", "w", encoding="utf-8")

for index, line in enumerate(file):
    line = line.strip("\n")
    data = line.split("\t")

    gene_id = data[2]
    gene_symbol = data[5]
    aliases = data[6]
    output_file.write(f"{gene_id}\t{gene_symbol}\t{aliases}\n")

output_file.close()