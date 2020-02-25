import argparse
from Bio import Entrez
import time
import json



## arg parse section ##
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--virus", help="type in a virus to search NCBI for", type=str)
parser.add_argument("-e", "--email", help="please include your email for NCBI", type=str)
args = parser.parse_args()

virus = args.virus
email = args.email

out_json = "%s_NCBI_results.json" % virus
out_fasta = "%s_NCBI_results.fasta" % virus
out_fasta_errors = "%s_NCBI_results_errors.fasta" % virus


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
    print(f"Gathering item {count+1} out of {len(acc_nums)}")


with open(out_json, "w") as j_out:
    json.dump(data, j_out, indent=2)



with open(out_json, "r") as in_j:
    in_data = json.load(in_j)

print("DONE GATHERING METADATA, WRITING {virus} INFORMATION TO FASTA FILE")

fasta = []
for item, value in in_data.items():
    temp_host = value["GBSeq_feature-table"][0]["GBFeature_quals"]
    host = []
    for i in temp_host:
        if i["GBQualifier_name"] == "host":
            host.append(i["GBQualifier_value"])
    if len(host) == 0:
        host = ["NHI"]

    date = value["GBSeq_create-date"]
    if date == None:
        date = "ND"


    header = "%s_%s_%s" %(item+"_"+virus, host[0], date)
    header = header.replace(".", "_")
    header = header.replace("-", "_")
    header = header.replace(" ", "_")
    header = header.replace(";","_")
    header = header.replace(",", "_")
    #print(header)

    seq = list(value["GBSeq_sequence"])
    if len(seq) >= 20000 and len(seq) <=50000:
        seq = "".join(seq)
    else:
        seq = "error"

    fasta.append((header, seq))

with open(out_fasta, "w") as out_f, open(out_fasta_errors, "w") as out_e:
    for line in fasta:
        if line[1] == "error":
            out_e.write(">{}\n{}\n".format(line[0],line[1]))
        else:
            out_f.write(">{}\n{}\n".format(line[0],line[1]))

print(f"ALL {virus} GENOMES WRITTEN TO FASTA FILE")
