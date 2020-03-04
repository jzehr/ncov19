import argparse
from Bio import Entrez
import time
import json
import csv
from collections import Counter


from helper import jde, fasta_writer

## arg parse section ##
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--virus", help="type in a virus to search NCBI for", type=str)
parser.add_argument("-e", "--email", help="please include your email for NCBI", type=str)
args = parser.parse_args()

virus = args.virus
email = args.email

out_json = "data/jsons/%s_NCBI_results.json" % virus
out_fasta = "data/fasta/%s_NCBI_results.fasta" % virus
out_fasta_errors = "data/fasta/%s_NCBI_results_errors.fasta" %  virus


Entrez.email = email

search_term = virus + " coronavirus complete genome"

print(f"SEARCHING FOR {virus} ON NCBI\n")

search = Entrez.esearch(db="nucleotide", retmax=10000, term=search_term, idtype="acc")
records = Entrez.read(search)

acc_nums = records["IdList"]
print(f"FOUND {len(acc_nums)} {virus} GENOMES\n")


print("SEARCHING FOR METADATA... BE PATIENT\n")

data = {}
for count, acc in enumerate(acc_nums):
    temp = {}
    handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="xml")
    record = Entrez.read(handle)[0]
    temp[acc] = record

    data.update(temp)
    time.sleep(1)
    print(f"Gathering genome {count+1} out of {len(acc_nums)}")

with open(out_json, "w") as j_out:
    json.dump(data, j_out, indent=2)

with open(out_json, "r") as in_j:
    in_data = json.load(in_j)


keys = list(in_data.keys())
info = {}
for key in keys:
    info.update(jde(key, in_data))


good_info = {}
error_info = {}

## write master csv here ##
# master row should contain all key information from good_info #
virus_data_csv = "data/csvs/%s_info.cvs" %  virus
with open(virus_data_csv, "w") as out_csv:
    spamwriter = csv.writer(out_csv, delimiter="\t")
    row_key = list(info.keys())[0]
    v = list(info[row_key].keys())
    row = ["ACC_NUM"] + v


    spamwriter.writerow(row)

    for key, value in info.items():
        if type(key) != tuple:

            temp = {}
            temp[key] = value
            spamwriter.writerow( [key] + list(list(temp.values())[0].values()))
            good_info.update(temp)

        else:
            temp = {}
            temp[key[1]] = value
            error_info.update(temp)

good_outs = "data/jsons/%s_good_info.json" % virus
with open(good_outs, "w") as j_out:
    json.dump(good_info, j_out, indent=2)

bad_outs = "data/jsons/%s_bad_info.json" % virus
with open(bad_outs, "w") as j_out:
    json.dump(error_info, j_out, indent=2)
print(f"DONE GATHERING METADATA, WRITING {virus} INFORMATION TO FASTA FILE")


with open(good_outs, "r") as j_in:
    meta_data = json.load(j_in)

## need to cylce through and grab all genes ##

keys = list(meta_data.keys())

genes = []
for key in keys:
    genes.append(meta_data[key]["genes"])

new_genes = list(Counter([i for j in genes for i in j]))

virus_gene_json = {virus : new_genes}
vgj = "data/jsons/%s_genes.json" % virus
with open(vgj, "w") as outer:
    json.dump(virus_gene_json, outer, indent=2)

for g in new_genes:
    ## write file here ##
    virus_gene_fasta = "data/fasta/%s_%s.fasta" % (virus, g)
    with open(virus_gene_fasta, "w") as out:
        for key in keys:
            list_genes = meta_data[key]["genes"]
            list_index = meta_data[key]["indices"]
            for pos, item in enumerate(list_genes):
                if item == g:
                    index = list_index[pos]
                    results = fasta_writer(meta_data, key, index, g)
                    header, seq, row = results[0], results[1], results[2]
                    out.write(">{}\n{}\n".format(header, seq))




