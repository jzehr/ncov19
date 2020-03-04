import json


SARS_gene_json = "data/jsons/SARS_genes.json"
MERS_gene_json = "data/jsons/MERS_genes.json"


with open(SARS_gene_json, "r") as sars_in, open(MERS_gene_json, "r") as mers_in:
  sars_data, mers_data = json.load(sars_in), json.load(mers_in)

sars_g_list, mers_g_list = list(sars_data.values())[0], list(mers_data.values())[0]
s, m = ["SARS"]*len(sars_g_list), ["MERS"]*len(mers_g_list)

sars, mers = [i + "_"+ j for i, j in zip(s, sars_g_list)], [i + "_"+ j for i, j in zip(m, mers_g_list)] 

total_files = sars + mers

#print(total_files)


### READY FOR PIPELINE ###

#HYPHY = "hyphy-develop/hyphy LIBPATH=hyphy-develop/res"
PRE = "hyphy-analyses/codon-msa/pre-msa.bf"
POST = "hyphy-analyses/codon-msa/post-msa.bf"
GARD = "hyphy-develop/res/TemplateBatchFiles/GARD.bf"
MEME = "hyphy-develop/res/TemplateBatchFiles/SelectionAnalyses/MEME.bf"
FEL_contrast = "hyphy-develop/res/TemplateBatchFiles/SelectionAnalyses/FEL-contrast.bf"


rule all:
  input:
    expand("data/fasta/{vir_seq}.fasta_protein_aligned.fas", vir_seq=total_files)


####################################################################
# This rule will read in a reg_prod_file and run it through pre-bf
# ~~ some proteins are not in frame... ~~ might want to pre-process 
# this data to get rid of these first 
####################################################################
rule rpf_pre:
  input:
    in_f = "data/fasta/{vir_seq}.fasta"
  output:
    out_prot = "data/fasta/{vir_seq}.fasta_protein.fas",
    out_nuc = "data/fasta/{vir_seq}.fasta_nuc.fas"
  shell:
   "hyphy-develop/hyphy {PRE} --input {input.in_f} --E 0.05"

####################################################################
# This rule will read protein fas from previous rule and
# align it with MAFFT
#
# add a separate dir to send these files
####################################################################
rule mafft_rpf:
  input:
    in_prot = rules.rpf_pre.output.out_prot
  output:
    out_prot = "data/fasta/{vir_seq}.fasta_protein_aligned.fas"
  shell:
    "mafft --quiet {input.in_prot} > {output.out_prot}"

####################################################################
# This rule will read in aligned PROTEIN file
# and run it through post-bf
####################################################################
rule rpf_post:
  input:
    in_prot = rules.mafft_rpf.output.out_prot,
    in_nuc = rules.rpf_pre.output.out_nuc
  output:
    out_f = "data/fasta/{vir_seq}.fasta_protein_aligned.fas.hyphy.fas"
  shell:
   "hyphy-develop/hyphy {POST} --protein-msa {input.in_prot} --nucleotide-sequences {input.in_nuc} --output {output.out_f} --compress No"

####################################################################
# This rule will read in the post-hyphy fasta 
# and run it through ~ GARD ~
####################################################################
rule rpf_GARD:
  input:
    in_f = rules.rpf_post.output.out_f
  output:
    out_j = str(rules.rpf_post.output.out_f) + ".GARD.json",
    out_nex = str(rules.rpf_post.output.out_f) + ".best-gard"
  shell:
   "hyphy-develop/hyphy {GARD} --alignment {input.in_f}"
