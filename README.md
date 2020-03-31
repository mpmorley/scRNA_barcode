# scRNA_barcode
```bash
usage: ProccessReads.py [-h] -1 READ1 -2 READ2 -s STAGGERLENGTH
                        [-g GFPPRIMERLENGTH] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -1 READ1, --read1 READ1
                        1st Read Pair one with Cellid+UMI
  -2 READ2, --read2 READ2
                        2nd Read Pair
  -s STAGGERLENGTH, --staggerLength STAGGERLENGTH
                        Stagger Length
  -g GFPPRIMERLENGTH, --gfpPrimerLength GFPPRIMERLENGTH
                        length of GFP primer
  -o OUTPUT, --output OUTPUT
                        Name of the output directory

python ProccessReads.py -1 AATCCAGC_1.fastq.gz -2 AATCCAGC_2.fastq.gz -s 4 -o Sample1

```

The next step is run using the R script ProcessBarcode.R

You'll need to have the starcode binary in your sytem path. 

This fiorst thing to do is configure the parameter sweep you will run with starcode. 
This done by modifing the samplerRun listat hte top of the file. 
```bash

samplerRun<- list()
samplerRun[['sub50_d8']] <- list(subSampleSize=50, d=8)
samplerRun[['sub40_d8']] <- list(subSampleSize=40, d=8)
samplerRun[['sub30_d8']] <- list(subSampleSize=30, d=8)
samplerRun[['sub30_d6']] <- list(subSampleSize=30, d=6)
```
Then run the script. The Whitelist.tsv file is the Whitelist cell barcodes, make sure it matches the chemistry version that was used for the run. 
```bash
Rscript --vanilla ProcessBarcode.R -s Sample1/uniqueShavedReads.txt -w Whitelist.tsv -o Sample

```




