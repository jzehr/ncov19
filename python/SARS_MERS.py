import argparse
from Bio import Entrez
import time
import json
from collections import Counter


from helper import jde, fasta_writer

## arg parse section ##
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--virus", help="type in a virus to search NCBI for", type=str)
parser.add_argument("-e", "--email", help="please include your email for NCBI", type=str)
args = parser.parse_args()

virus = args.virus
email = args.email

out_json = "%s/%s_NCBI_results.json" % (virus, virus)
out_fasta = "%s/%s_NCBI_results.fasta" % (virus, virus)
out_fasta_errors = "%s/%s_NCBI_results_errors.fasta" % (virus, virus)


Entrez.email = email

search_term = virus + " coronavirus complete genome"

print(f"SEARCHING FOR {virus} ON NCBI\n")

'''
search = Entrez.esearch(db="nucleotide", retmax=10000, term=search_term, idtype="acc")
records = Entrez.read(search)

acc_nums = records["IdList"]
print(f"FOUND {len(acc_nums)} {virus} GENOMES\n")


print("SEARCHING FOR METADATA... BE PATIENT\n")
'''

'''
data = {}
for count, acc in enumerate(acc_nums):
    temp = {}
    handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="xml")
    record = Entrez.read(handle)[0]
    temp[acc] = record

    data.update(temp)
    time.sleep(1)
    print(f"Gathering item {count+1} out of {len(acc_nums)}")

with open(out_json, "w") as j_out:
    json.dump(data, j_out, indent=2)
'''

with open(out_json, "r") as in_j:
    in_data = json.load(in_j)


keys = list(in_data.keys())
info = {}
for key in keys:
    info.update(jde(key, in_data))

good_info = {}
error_info = {}
for key, value in info.items():
    if type(key) != tuple:
        temp = {}
        temp[key] = value
        good_info.update(temp)
    else:
        temp = {}
        temp[key[1]] = value
        error_info.update(temp)

good_outs = "%s/%s_good_info.json" % (virus, virus)
with open(good_outs, "w") as j_out:
    json.dump(good_info, j_out, indent=2)

bad_outs = "%s/%s_bad_info.json" % (virus, virus)
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

for g in new_genes:
    ## write file here ##
    virus_gene = "%s/%s_%s.fasta" % (virus, virus, g)
    with open(virus_gene, "w") as out:
        for key in keys:
            list_genes = meta_data[key]["genes"]
            list_index = meta_data[key]["indices"]
            for pos, item in enumerate(list_genes):
                if item == g:
                    index = list_index[pos]
                    results = fasta_writer(meta_data, key, index, g)
                    header, seq, row = results[0], results[1], results[2]
                    out.write(">{}\n{}\n".format(header, seq))
## then open a file for each gene to write seqs to ##



'''
with open(out_fasta, "w") as out_f, open(out_fasta_errors, "w") as out_e:
    for line in fasta:
        if line[1] == "error":
            out_e.write(">{}\n{}\n".format(line[0],line[1]))
        else:
            out_f.write(">{}\n{}\n".format(line[0],line[1]))

print(f"ALL {virus} GENOMES WRITTEN TO FASTA FILE")
'''
