import json
from nested_lookup import nested_lookup
from collections import Counter
import re
import math



def name_fixer(string):
    new = []
    bad_chars = [" ",",",".","/","-","(",")",":","'", ";"]
    def checker(char):
        if char in bad_chars:
            return "_"
        else:
            return char
    for char in list(string):
        new.append(checker(char))

    temp = "".join(new)
    name = temp.replace("__","_")
    return name

def checker(lst):
        temp = lst
        if len(temp) == 0:
            temp.append("no_value")
        return temp

def jde(key, dictionary):

    master = {}
    data = {}

    error = {}

    ## accesion numbers ##
    accession = name_fixer(str(key))
    value = dictionary[key]

    ## setting return dict ##
    master[accession] = data

    ## ~~ name of virus ~~ ##
    virus_name = name_fixer(checker(value["GBSeq_definition"]))


    ## ~~ create data for seq ~~ ##
    create_date = name_fixer(checker(value["GBSeq_create-date"]))

    ## ~~ paper ~~ ##
    paper_title = name_fixer(checker(value["GBSeq_references"][0]["GBReference_title"]))

    ## doi for paper if available ##
    doi = []
    temp_doi = value["GBSeq_references"][0]#["GBReference_xref"]
    if "GBReference_xref" in list(temp_doi.keys()):
        doi.append(temp_doi["GBReference_xref"][0]["GBXref_id"])
    else:
        doi.append("no_doi")

    ## ~~ to return ~~ ##
    doi = name_fixer(checker(doi[0]))

    other_quals = value["GBSeq_feature-table"][0]["GBFeature_quals"]
    host = []
    isolate = []
    country = []
    collection_date = []
    for d in other_quals:
        info = list(d.values())
        if info[0] == "host":
            host.append(info[1])
        elif info[0] == "isolate":
            isolate.append(info[1])
        elif info[0] == "country":
            country.append(info[1])
        elif info[0] == "collection_date":
            collection_date.append(info[1])
        else:
            continue

    ## ~~ to return ~~ ##
    host = name_fixer(checker(host))
    isolate = name_fixer(checker(isolate))
    country = name_fixer(checker(country))
    collection_date = name_fixer(checker(collection_date))

    feat_tab = nested_lookup("GBSeq_feature-table", value)
    ## ~~ to return ~~ ##
    gene_names = []
    locations = []
    trans_seq = []
    for feat in feat_tab:
        for i in feat:

            ## getting location of seq ##
            if i["GBFeature_key"] == "gene":
                locations.append(i["GBFeature_location"])

            ## getting name of gene
            elif i["GBFeature_key"] == "CDS":
                interval = i["GBFeature_intervals"]
                quals = i["GBFeature_quals"]

                for j in quals:
                    ## gene names ##
                    if j["GBQualifier_name"] == "gene":
                        gene_names.append(j["GBQualifier_value"])

                    ## translated seq ##
                    elif j["GBQualifier_name"] == "translation":
                        trans_seq.append(j["GBQualifier_value"])

    # sanity check to make sure genes, incides, and trans seqs match in length #
    #print(len(gene_names), len(locations), len(trans_seq))

    ## whole sequence ##
    seq = value["GBSeq_sequence"]

    new_gene_names = []
    for pos, gene in enumerate(gene_names):
        q = gene.upper() + "_" + str(pos)
        if q == "ORF1AB_0":
            new_gene_names.append(gene.upper())
        elif q == "ORF1AB_1":
            new_gene_names.append("ORF1A")
        else:
            new_gene_names.append(gene.upper())

    data["virus_name"] = virus_name
    data["host"] = host
    data["isolate"] = name_fixer(checker(isolate))
    data["country"] = country
    data["collection_date"] = collection_date
    data["paper_title"] = paper_title
    data["doi"] = doi
    data["create_date"] = create_date


    ## this could be fragile, but it ensures that each gene has a location ##
    if len(gene_names) == len(locations):
        data["genes"] = new_gene_names
        data["indices"] = locations
        data["trans_seq"] = trans_seq
        data["nuc_seq"] = seq
        return master

    else:
        error[("prodCDS_length", name_fixer(accession))] = data
        return error

def fasta_writer(meta, key, index, gene):
    acc_num = key
    nuc_seq = meta[key]["nuc_seq"]
    orf1_s = int(re.findall(r'\d+', index.split("..")[0])[0]) - 1
    orf1_e = int(re.findall(r'\d+', index.split("..")[1])[0])
    date = []
    if meta[key]["collection_date"][0] == "no_value":
        date.append(("create", meta[key]["create_date"]))
    else:
        date.append(("collection_date", meta[key]["collection_date"]))

    country = name_fixer(meta[key]["country"])
    host = name_fixer(meta[key]["host"])
    isolate = name_fixer(meta[key]["isolate"])
    date_type = date[0][0]
    date_date = name_fixer(date[0][1])
    nuc_seq = meta[key]["nuc_seq"]

    #print(f"host {host} | country {country} | isolate {isolate} | date_type {date_type} | date_date {date_date}")

    seq = nuc_seq[orf1_s:orf1_e]
    header = "%s_%s_%s_%s_%s_%s_%s" % (name_fixer(acc_num), gene, date_type, date_date, name_fixer(country), name_fixer(host), name_fixer(isolate))
    header = header.replace("__","_")
    row = (name_fixer(acc_num), name_fixer(gene), date_type, date_type, name_fixer(country), name_fixer(host), name_fixer(isolate))
    return header, seq, row


def distance(location_1, location_2):
    lat1, lon1 = location_1
    lat2, lon2 = location_2
    radius = 6371

    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) * math.sin(dlat / 2) +
            math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
	    math.sin(dlon / 2) * math.sin(dlon / 2))
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    d = radius * c
    return d

