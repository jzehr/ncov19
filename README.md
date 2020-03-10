REQUIREMENTS:
python 3 or higher
Anaconda3

## Step 1: Get hyphy-analyses
```
git clone https://github.com/veg/hyphy-analyses.git
git clone https://github.com/veg/hyphy.git hyphy-develop
```
### HYPHY will be installed using conda (hyphy v. 2.5.5)

## Step 2: start conda
```
cd ..

conda env create -f environment.yml

conda activate viran
```

### Might have to:
```
unset PYTHONPATH

pip install nested-lookup
```

## Step 3: gather information (only need to do this once!)
```
python python/SARS_MERS.py -e [EMAIL] -v [VIRUS]

python python/SARS_MERS.py -e gmail@gmail.com -v SARS
```


## Step 4: Run the pipeline with snakemake
### This is where things get fun
```
bpsh 3 snakemake data/fasta/MERS_SARS_COVID19_S.fasta_protein_aligned.fas.hyphy.fas.best-gard.nex.ALL.json -j 200

(or one particular analysis)

bpsh 3 snakemake data/fasta/MERS_SARS_COVID19_S.fasta_protein_aligned.fas.hyphy.fas.best-gard.nex.MEME.json -j 200
```

-- you can use a ```-j``` flag to denote the number of cores to run on.

If you are on the server:

-- try to use ```bpsh 2``` to send the process somewhere besides the head node

