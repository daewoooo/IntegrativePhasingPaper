This repository contains the pipeline used to generate the results reported in the paper ["Dense And Accurate Whole-Chromosome Haplotyping Of Individual Genomes"](https://doi.org/10.1101/126136). Please also refer to that paper for a description of the method.

To install all dependencies, we recommend [Conda](http://anaconda.org). If you don't yet use conda, it is usually easiest to just install [miniconda](https://conda.io/miniconda.html),  see [instructions](https://conda.io/docs/install/quick.html).

Then, create and activate an environment that has all required packages installed (where ``environment.yml`` is the file we provide here).

```
conda env create --name integrative-phasing -f environment.yml
source activate integrative-phasing
snakemake --resources download=8 -j 40
```

This command will run the full pipeline using 40 threads and at most 8 downloads at a time.

