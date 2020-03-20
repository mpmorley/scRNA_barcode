# scRNA_barcode

usage: ProccessReads.py [-h] -1 READ1 -2 READ2 -s STAGGERLENGTH
                        [-g GFPPRIMERLENGTH] [-o OUTPUT]
optional arguments:   -h, --help            show this help message and exit   -1 READ1, --read1 READ1
                        1st Read Pair one with Cellid+UMI   -2 READ2, --read2 READ2
                        2nd Read Pair   -s STAGGERLENGTH, --staggerLength STAGGERLENGTH
                        Stagger Length   -g GFPPRIMERLENGTH, --gfpPrimerLength GFPPRIMERLENGTH
                        length of GFP primer   -o OUTPUT, --output OUTPUT
                        Name of the output directory



> python ProccessReads.py -1 AATCCAGC_1.fastq.gz -2 AATCCAGC_2.fastq.gz
> -s 4 -o Sample1


