# Install

After cloning the github repository, you will need to install the conda environment specified in the .yaml file.

```conda env create -f blast.yaml
conda activate blast
```

# Usage

To use the package, you will need to have your sequence in a fasta file, and your database in another. 
Run from the root repository

```python3 src/blast.py data/name_of_sequence_file.fasta data/name_of_db_file.fasta```

You can test out the package with the sample data : 

```python3 src/blast.py data/P06007.fasta data/photosysII_chlamydomonas.fasta > results/chlamydomonasdb.txt```

# Search parameters

Currently, the search parameters can be changed in the blast.py script in the 7 to 14 first lines of the code

# Output

Currently the results will be in the default output. You can thus redirect it to a file by using the ">" operator as shown in the following example : 

```python3 src/blast.py data/name_of_sequence_file.fasta data/name_of_db_file.fasta > results/result_file.txt```
