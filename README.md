This repository contains the pipeline used to generate the results reported in the paper ["Dense And Accurate Whole-Chromosome Haplotyping Of Individual Genomes"](https://doi.org/10.1101/126136). Please also refer to that paper for a description of the method.

To install all dependencies, we recommend [Conda](http://anaconda.org). If you don't yet use conda, it is usually easiest to just install [miniconda](https://conda.io/miniconda.html),  see [instructions](https://conda.io/docs/install/quick.html).

Then, create and activate an environment that has all required packages installed (where ``environment.yml`` is the file we provide here).

```
conda env create --name integrative-phasing -f environment.yml
source activate integrative-phasing
snakemake --resources download=8 -j 40
```

This command will run the full pipeline using 40 threads and at most 8 downloads at a time.

Note that we ran a lot of different combinations of parameters for our paper. Please edit the first few lines of the Snakefile according to how many experiments you want to run. Here's the corresponding snippet:

```
#coverage = [2,3,4,5,10,15,25,30,'all']
coverage = ['all']
#chromosome = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
chromosome = [1]
#strandseqcoverage = [5,10,20,40,60,80,100,120,134]
strandseqcoverage = [134]
#trials = [1,2,3,4,5]
trials = [3]
```
