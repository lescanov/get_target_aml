### Purpose:
* Retrieve raw RNA and miRNA-seq counts from TARGET-AML cohort, specifically filtering for:
    * Bone marrow samples
    * First replicate, if many replicates exist
* Format clinical data:
    * Clean columns for easier access
    * Subcategorize subtypes, spefifically the pre-defined "Other" subgroup into t(6;9), t(3;5)(q25;q34), del(5q), del(7q), del(9q), monosomy 5 or 7, trisomy 8 or 21, monosomy Y or X

This workflow uses TCGAbiolinks to retrieve raw data from the GDC. Clinical data was retrieved from cbioportal.

### How to use:
You must have snakemake installed. Fork the repository and use the following code:

If using a non-ARM based machine:
```zsh
snakemake --use-conda --cores 1
```

If using ARM based machine:
```zsh
CONDA_SUBDIR=osx-64 snakemake --use-conda --cores 1
```

** Depending on your machine specs, this workflow may crash **  
Retrieving raw counts for TARGET will use lots of ram. May need to modify your Renv to allow R to utilize more resources. If there is a problem, you will have to modify the amount of RAM specified to use in the snakemake rules pertaining to the step in question.
