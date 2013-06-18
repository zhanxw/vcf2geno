vcf2geno
=========

Convert VCF files to genotype file and site file.

Input File
----------

Vcf2geno takes VCF files. They can be in plain text or GZIP/BGZIP compressed formats.

Outputs
-------
Vcf2geno generates two sets of files: prefix.geno and prefix.site where ''prefix'' is the given parameter to ''--out''.

A .geno file is shown below:

    X1  X1  0   -9  0   0   0   2   2   2   2   2   2
    X2  X2  -9  -9  0   0   0   2   2   2   2   2   2
    X3  X3  0   -9  0   0   0   2   2   2   2   2   2
    X4  X4  0   -9  0   0   0   2   2   2   2   2   2
    X5  X5  0   -9  0   0   0   2   2   1   2   2   2
    X6  X6  0   -9  0   0   0   2   2   2   2   2   -9

The first and second columns are sample IDs copied from the header of VCF files.
From column 3 till the last the column, they are individual level genotype converted from VCF files.

A .site file is shown below:

    CHROM   POS ID  REF ALT
    1   10  1:10    A   T
    1   20  1:20    G   C
    1   30  1:30    C   A
    1   40  1:40    A   C
    1   10000   1:10000 G   C
    1   20000   1:20000 T   A
    4   5000    4:5000  A   T
    4   6000    4:6000  C   T
    X   800 X:800   A   C
    X   900 X:900   A   T
    X   1000    X:1000  T   G

The content of .site file begins with a header line, and thus the content part from the second line is chromosome, positions, reference alleles and alternative alleles.

Options
-------

Vcf2geno provides sample selection options and range selection options.

There are four options to include/exclude samples:

    --peopleIncludeID, --peopleExcludeID: specify which samples are included/excluded in conversion 

e.g. --peopleIncludeID X1,X2,X3 will convert only 3 people during conversion if input VCF file contains these three samples.

    --peopleIncludeFile, --peopleExcludeFile: specify a file to include/exclude samples. 
    
Each line of the file should be a sample ID.

There are two options to specify regions. You can convert part of the VCF file using this option, however, your input file must be indexed by TABIX.

    --rangeList: this options enable you to specify a range by hand. 
e.g. --rangeList 1:100-200. Note your chromosome name in the command line should be consistent to the content of the VCF file (e.g. both do not have 'chr' prefix).

    --rangeFile: this option specifies range using an external file. 
Each line of the file should specify a range, e.g. '1:100-200', or alternatively three columns '1 100 200'. 

    --updateID: this option specifies an external file to update the sample ID in the input VCF file.
This file should have three columns. The first columns are sample ID used in the input VCF file. The second and third columns are the family ID and individual ID. 
vcf2geno will update the IDs in the VCF to the family ID and individual ID in the outputted genotype file.
 
Example
-------

Under the "exampleVCF" folder, you can file example.vcf.gz. This is an indexed VCF file.
Basic usage of extracting all samples across all regions:

    ../vcf2geno --inVcf example.vcf.gz --out test

Convert sample X1 from range 1:20-30:

    ../vcf2geno --inVcf example.vcf.gz --rangeList 1:20-30 --peopleIncludeID X1 --out test

Extract all samples and update sample IDs using a separate file ''test.updateId'':

    ../vcf2geno --inVcf example.vcf.gz --out test --updateID test.updateId 

 
Contact
-------
Questions or comments should be sent to Xiaowei Zhan
([zhanxw@umich.edu](mailto:zhanxw@umich.edu "mailto:zhanxw@umich.edu"))
or Chaolong Wang
([chaolong@umich.edu](mailto:chaolong@umich.edu "mailto:chaolong@umich.edu"))
or Goncalo Abecasis
([goncalo@umich.edu](mailto:goncalo@umich.edu "mailto:goncalo@umich.edu"))


