# CAS-523-Proj-2

## Genome

### Wuhan-Hu-1
We are using the complete genome of the Sars-COV-2 virus strain Wuhan-Hu-1
found [National Library of Medicine](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512)
for the first part of the project. The text version of the genome can be
found in `original/`.

When reading this genome, we can break it into the 5'UTR region (sites 1-265) 
and the actual gene (sites 266-21555) according to the FASTA file linked in the 
above reference. 

The `original/fasta_reader.py` can be used to read the `wuhan-hu-1.txt` file 
by calling its `read_wuhan_1` function with the filepath to `wuhan-hu-1.txt`.
This function will return 3 strings: `hu1_full_genome, utr5, hu1_gene`. The 
full genome is just that. `utr5` is the 5'UTR region (first 265 nucleotides),
and `hu1_gene` is a string comprised of the nuclotides at sites 266-21555. 
It has an arg parser and can be called in the terminal:

```shell
python fasta_reader.py [-filepath=<file/you/are/reading> -opt=[which file 
you are reading]
```

The only option for `-opt` at this stage is `0` for ``read_wuhan_1`. If we 
use a different FASTA for part 2, we can update the reader to portion it 
appropriately and add its option. 
