import json
import os

SARS_gene_json = "data/jsons/SARS_genes.json"
MERS_gene_json = "data/jsons/MERS_genes.json"


with open(SARS_gene_json, "r") as sars_in, open(MERS_gene_json, "r") as mers_in:
  sars_data, mers_data = json.load(sars_in), json.load(mers_in)

sars_g_list, mers_g_list = list(sars_data.values())[0], list(mers_data.values())[0]
s, m = ["SARS"]*len(sars_g_list), ["MERS"]*len(mers_g_list)

sars, mers = [i + "_"+ j for i, j in zip(s, sars_g_list)], [i + "_"+ j for i, j in zip(m, mers_g_list)] 

total_files = sars + mers

total_files.remove("SARS_ORF1AB")
total_files.remove("SARS_PP1AB")
total_files.remove("MERS_PP1AB")
total_files.remove("MERS_ORF1AB")

#print(total_files)


### READY FOR PIPELINE ###

PRE = "hyphy-analyses/codon-msa/pre-msa.bf"
POST = "hyphy-analyses/codon-msa/post-msa.bf"

GARD = "hyphy-develop/res/TemplateBatchFiles/GARD.bf"
MEME = "hyphy-develop/res/TemplateBatchFiles/SelectionAnalyses/MEME.bf"
#FEL_contrast = "hyphy-develop/res/TemplateBatchFiles/SelectionAnalyses/FEL-contrast.bf"
FEL = "hyphy-develop/res/TemplateBatchFiles/SelectionAnalyses/FEL.bf"
FUBAR = "hyphy-develop/res/TemplateBatchFiles/SelectionAnalyses/FUBAR.bf"
SLAC = "hyphy-develop/res/TemplateBatchFiles/SelectionAnalyses/SLAC.bf"
FMM = "hyphy-analyses/FitMultiModel/FitMultiModel.bf"
BUSTED = "hyphy-develop/res/TemplateBatchFiles/SelectionAnalyses/BUSTED.bf"

rule all:
  input:
    expand("data/fasta/{vir_seq}.fasta_protein_aligned.fas.hyphy.fas.GARD.json", vir_seq=total_files)


####################################################################
# This rule will read in a reg_prod_file and run it through pre-bf
# ~~ some proteins are not in frame... ~~ might want to pre-process 
# this data to get rid of these first 
####################################################################
rule pre_mafft:
  input:
    in_f = "data/fasta/{vir_seq}.fasta"
  output:
    out_prot = "data/fasta/{vir_seq}.fasta_protein.fas",
    out_nuc = "data/fasta/{vir_seq}.fasta_nuc.fas"
  shell:
   "hyphy {PRE} --input {input.in_f} --E 0.001"

####################################################################
# This rule will read protein fas from previous rule and
# align it with MAFFT
#
# add a separate dir to send these files
####################################################################
rule mafft_align:
  input:
    in_prot = rules.pre_mafft.output.out_prot
  output:
    out_prot = "data/fasta/{vir_seq}.fasta_protein_aligned.fas"
  shell:
    "mafft --quiet {input.in_prot} > {output.out_prot} 2> mafft_errors.log"

####################################################################
# This rule will read in aligned PROTEIN file
# and run it through post-bf
####################################################################
rule post_mafft:
  input:
    in_prot = rules.mafft_align.output.out_prot,
    in_nuc = rules.pre_mafft.output.out_nuc
  output:
    out_f = "data/fasta/{vir_seq}.fasta_protein_aligned.fas.hyphy.fas"
  shell:
   "hyphy {POST} --protein-msa {input.in_prot} --nucleotide-sequences {input.in_nuc} --output {output.out_f} --compress No"

####################################################################
# This rule will read in the post-hyphy fasta 
# and run it through ~ GARD ~
####################################################################
rule hyphy_GARD:
  input:
    in_f = rules.post_mafft.output.out_f
  output:
    out_j = str(rules.post_mafft.output.out_f) + ".GARD.json",
    out_nex = str(rules.post_mafft.output.out_f) + ".best-gard"
  shell:
    "hyphy {GARD} --alignment {input.in_f}"
 
########################################################################
# This rule will read in the post-GARD .best-gard file 
# if file is empty it will build a tree, otherwise it will use the file
########################################################################
rule tree_maker:
  input:
    in_nex = rules.hyphy_GARD.output.out_nex,
    in_aligned = rules.post_mafft.output.out_f
  output:
    out_nex = str(rules.hyphy_GARD.output.out_nex) + ".nex", 
    #out_tmp = str(rules.hyphy_GARD.output.out_nex) + ".tmp" 
  run:
    if not os.stat(input.in_nex).st_size == 0:
      shell("mv {input.in_nex} {output.out_nex}")
    #else:
    #  # build a tree then make a nexus
    #  shell("fasttree < {input.in_aligned} > {output.out_tmp}")
    #  shell("cat {input.in_aligned} {output.out_tmp} > {output.out_nex}")

#####################################################################
# This rule will read in the output of tree_maker 
# and run it through ~ MEME ~
####################################################################
rule hyphy_MEME:
  input:
    in_nex = rules.tree_maker.output.out_nex
  output:
    out_j = str(rules.tree_maker.output.out_nex) + ".MEME.json"
  #run:
  #  import pdb;pdb.set_trace()
  shell:
    "hyphy {MEME} --alignment {input.in_nex}" 

###################################################################
# This rule will read in the output of tree_maker 
# and run it through ~ FEL ~
####################################################################
rule hyphy_FEL:
  input:
    in_nex = rules.tree_maker.output.out_nex
  output:
    out_j = str(rules.tree_maker.output.out_nex) + ".FEL.json"
  shell:
    "hyphy {FEL} --alignment {input.in_nex}" 

####################################################################
# This rule will read in the output of tree_maker 
# and run it through ~ SLAC ~
####################################################################
rule hyphy_SLAC:
  input:
    in_nex = rules.tree_maker.output.out_nex
  output:
    out_j = str(rules.tree_maker.output.out_nex) + ".SLAC.json"
  shell:
    "hyphy {SLAC} --alignment {input.in_nex}" 

####################################################################
# This rule will read in the output of tree_maker 
# and run it through ~ BUSTED ~
####################################################################
rule hyphy_BUSTED:
  input:
    in_nex = rules.tree_maker.output.out_nex
  output:
    out_j = str(rules.tree_maker.output.out_nex) + ".BUSTED.json"
  shell:
    "hyphy {BUSTED} --alignment {input.in_nex}" 

####################################################################
# This rule will read in the output of tree_maker 
# and run it through ~ FUBAR ~
####################################################################
rule hyphy_FUBAR:
  input:
    in_nex = rules.tree_maker.output.out_nex
  output:
    out_j = str(rules.tree_maker.output.out_nex) + ".FUBAR.json"
  shell:
    "hyphy {FUBAR} --alignment {input.in_nex}" 

####################################################################
# This rule will read in the output of tree_maker 
# and run it through ~ FMM ~
####################################################################
rule hyphy_FMM:
  input:
    in_nex = rules.tree_maker.output.out_nex
  output:
    out_j = str(rules.tree_maker.output.out_nex) + ".FITTER.json"
  shell:
    "hyphy {FMM} --alignment {input.in_nex}" 

#################################################
# This rule will read in ALL selection analyses 
# and output master selection json
#################################################
rule all_selection:
  input:
    in_GARD = rules.hyphy_GARD.output.out_j,
    in_FMM = rules.hyphy_FMM.output.out_j,
    in_MEME = rules.hyphy_MEME.output.out_j,
    in_FUBAR = rules.hyphy_FUBAR.output.out_j,
    in_FEL = rules.hyphy_FEL.output.out_j,
    in_BUSTED = rules.hyphy_BUSTED.output.out_j,
    in_SLAC = rules.hyphy_SLAC.output.out_j
  output:
    out_j = str(rules.tree_maker.output.out_nex) + ".ALL.json"
  shell:
    "cat {input.in_GARD} {input.in_FMM} {input.in_MEME} {input.in_FUBAR} {input.in_FEL} {input.in_BUSTED} {input.in_SLAC} > {output.out_j}" 



