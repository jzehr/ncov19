REQUIREMENTS:
python 3 or higher
Anaconda3

## Step 1: Get hyphy-analyses
```
git clone https://github.com/veg/hyphy-analyses.git
git clone https://github.com/veg/hyphy.git hyphy-develop
cd hyphy-develop
git checkout develop
cmake ./
make -j MP
```

## Step 2: start conda
```
cd ..
conda env create -f environment.yml

conda activate viran
```

### Might have to:
```
unset PYTHONPATH

pip3 install nested-lookup
```

## Step 3: gather information
```
python3 SARS_MERS.py -e [EMAIL] -v [VIRUS]

python3 SARS_MERS.py -e gmail@gmail.com -v SARS
```


## Step 4: Run the pipeline with snakemake
```
snakemake data/fasta/SARS_S.fasta_protein_aligned.fas.hyphy.fas.GARD.json
```
-- you can use a ```-j``` flag to denote the number of cores to run on.

-- try to use ```bpsh 2``` to send the process somewhere besides the head node

