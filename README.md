# IMCA

## Improving Mapping with Contigs Assembly

When mapping reads from any NGS experiment
 (like ChIP-seq, DNase-seq etc.)
 there is always some loss of information due to the differences
 between the reference genome and the actual genome being sequenced.
 IMCA attempts to reduce this loss.
 In short, first contigs are assembled from sequenced reads -
 they should correspond to the actual sequence of the genome.
 Next, the contigs are aligned to the reference genome,
 to find out where those sequences come from.
 Then reads are mapped to them.
 If there are any reads that didn't map to the reference
 due to differences in sequence,
 it might be that they will map to some contig
 that will map to the reference;
 then we can transfer coordinates of this read's mapping
 to the reference.
 *(picture here would be nice)*

If you encounter any problems running IMCA
 or have any comments or suggestions,
 please let me know at a dot macioszek at mimuw.edu.pl.

### Instalation and requirements

No installation needed.
At least not for these scripts.
You might need to add executive rights using `chmod +x` command to the python scripts.
`mapped_segment.py` should be in the same directory as `IMCA.py`,
because `IMCA.py` imports it.
Also you might need to install some software needed to complete additional steps.
Specifically,
you need some mapper for your reads, mapper for contigs, and assembler to create contigs.
The script was tested using minimap2 and velvet,
but should work with output of any other software,
as long as formats are correct:
the mappings should be in .bam / .sam format
(preferalby with CIGAR strings assigned to every mapping),
and contigs in fasta format.

If you provide your files in bam/sam format, pysam package is needed.
 If you provide them in fasta/fastq format and you want IMCA to map them,
 mappy package is needed.

### Workflow

1. Map reads to the reference.
2. Construct contigs from the reads.
3. Potentially, filter contigs.
4. Map (align) contigs to the reference.
5. Map reads to the mapped contigs.
6. Transfer the mapping of the reads from contigs to reference.

Steps from 1 to 5 can be done with any software of your choice.
Either step 6 alone or steps 1, 4-6 can be done using `IMCA.py` script.
Below I present example workflow.

**1. Map reads to the reference.**

Using [minimap2](https://github.com/lh3/minimap2):

`minimap2 -a -x sr reference.fa reads.fq > reads2reference.sam 2> mapping.err`

You can also use any other mapper, like bowtie2, BWA, Tophat etc.
 You can also leave mapping to IMCA (or more specifically, to mappy - minimap2 wrapper).

**2. Construct contigs from the reads.**

Using [velvet](https://www.ebi.ac.uk/~zerbino/velvet/):

`velveth velvet_output 27 -fastq reads.fq`

`velvetg velvet_output`

This will generate directory `velvet_output` with several files,
including large Sequences, Roadmap, Graph and LastGraph, that you can remove
once `velvetg` ends running.
We are only going to need contigs.fa and potentially stats.txt for filtering.
`27` refers to the size of the k-mers constructed by velvet
and of course can be changed according to your needs.
Please refer to [velvet's documantation](https://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf)
for more details.

You can use any other assembler.
The important thing is to get contigs in fasta format in this step.

**3. Potentially, filter contigs.**

You might choose to filter out contigs that are too short,
or have too small coverage.
I encourage to look at the obtained contigs and decide what seems best
to obtain better results.

If you used velvet as assembler, you may use `filter_velvet_contigs.py` script.
It doesn't work with other assemblers, because it uses
velvet's naming convention.
It can filter by contig's length and coverage.
For example, run:

`./filter_velvet_contigs.py -i velvet_output/contigs.fa -o velvet_output/filtered_contigs.fa -l 200 -c 2`

to create file `filtered_contigs.fa`, that will contain
these contigs from `velvet/contigs.fa` that were at least 200 bp long
and had coverage higher or equal to 2.

**4. Map (align) contigs to the reference.**

Using minimap2:

`minimap2 -a -x asm10 reference.fa velvet_output/contigs.fa > contigs2reference.sam 2> mapping_contigs.err`

You can also use any other mapper capable of mapping large contigs
or tool designed specifically for aligning genomes,
like MUMmer, BLAST etc.

Important note: the output of this step must be in sam or bam format.
If you use a software that generates anything else,
you need to convert it to sam / bam format before last step.
(For MUMmer, check out this [link](https://www.biostars.org/p/185384/).)
For better accuracy, CIGAR strings should be available for every alignment.
See documentation for details.
*(That is, see it when I post it. I didn't write it yet.)*

**5. Map reads to the mapped contigs.**

This step is simmilar to the step 1; we just use potentially filtered contigs instead of reference.
For example, using minimap2:

`minimap2 -a -x sr velvet_output/contigs.fa reads.fq > reads2contigs.sam 2> mapping2contigs.err`

But you can use Bowtie2, BWA, Tophat etc.

**6. Transfer the mapping of the reads from contigs to reference.**

Example command:

`./IMCA.py -r reads2reference.sam -c contigs2reference.sam --reads2contigs reads2contigs.sam -o merged_mapping.sam`

This will generate file `merged_mapping.sam`
 with all the reads that were mapped to the reference
 and reads mapped to the contigs which mapped to the reference.
 By default, the reads mapped to the reference stay unchanged,
 but you may choose to change their coordinates when possible by specifying `--mapped2ref transfer` argument.
 See documentation for more details.

Another example command:

`./IMCA.py -r reads2reference.sam -c contigs2reference.sam --reads2contigs reads2contigs.sam -o merged_mapping.sam --mapped2ref transfer --keep-unmapped-contigs`

This command will also generate `merged_mapping.sam` file,
 but this time:

1) all the reads that were mapped to reference and contigs
 will be transferred via contigs.
 That won't change the number of mapped reads, but can change their coordinates
 and can change the number of mappings;
 e.g., read A can have three mappings to reference (one primary, two secondary or supplementary)
 but only one mapping to contig, and the contig has one mapping to the reference;
 in the previous run there would be three mappings in the output,
 while in this one there would be one.

2) all the reads that were mapped to contigs that didn't map to reference
 will be outputted as mapped and the contigs will be added to header section of output,
 specifically to the section that describes references / contigs / chromosomes present in the SAM file.
 In the previous run such reads would be considered unmapped,
 so this argument can change the number of mapped reads.

### Shorter workflow

We can use IMCA to perform 1, 4-6 steps for us, using mappy.
 Even though it's a minimap2 wrapper,
 it gives slightly different results;
 also IMCA won't leave any .sam or .bam files except for the output,
 so you won't be able to compare how many reads were mapped where and so on.
 So, unless it's super incovenient for you to map everything by yourself,
 I recommend the above workflow.
 However if it *is* super inconvenient,
 you can use the workflow below.
 Steps 1 i 2 are the same as 2 and 3 from above, go there for more details.

1. Construct contigs from the reads.
2. Potentially, filter contigs.
3. Map everything and transfer the mapping of the reads from contigs to reference, all in one step.

**1. Construct contigs from the reads.**

Using [velvet](https://www.ebi.ac.uk/~zerbino/velvet/):

`velveth velvet_output 27 -fastq reads.fq`

`velvetg velvet_output`

You can use any other assembler.
The important thing is to get contigs in fasta format in this step.

**2. Potentially, filter contigs.**

If you used velvet as assembler, you may use `filter_velvet_contigs.py` script.
For example, run:

`./filter_velvet_contigs.py -i velvet_output/contigs.fa -o velvet_output/filtered_contigs.fa -l 200 -c 2`

to create file `filtered_contigs.fa`, that will contain
these contigs from `velvet/contigs.fa` that were at least 200 bp long
and had coverage higher or equal to 2.

**3. Map everything and transfer the mapping of the reads from contigs to reference.**

`./IMCA.py -r reads.fq -c velvet_output/filtered_contigs.fa --reference ref.fa -o merged_mapping.sam`

This will generate file `merged_mapping.sam`
 with all the reads that were mapped to the reference
 and reads mapped to the contigs which mapped to the reference.


### Parameters

To see all the arguments available, run the script with `-h` option:

```
usage: IMCA.py [-h] [-r READS] [-c CONTIGS] [--reference REFERENCE]
               [--reads2contigs READS2CONTIGS] [-o OUTPUT]
               [--mapped2ref MAPPED2REF] [-k]

You should always provide reads and contigs. If they are already mapped and in
sam/bam format, you should also provide reads2contigs. If they are in
fasta/fastq format and you want me to map them, you should also provide
reference.

optional arguments:
  -h, --help            show this help message and exit
  -r READS, --reads READS
                        reads; either bam/sam file with reads mapped to
                        reference or fasta/fastq file to be mapped to
                        reference with minimap.
  -c CONTIGS, --contigs CONTIGS
                        contigs; either bam/sam file with contigs mapped to
                        reference or fasta/fastq file to be mapped to
                        reference with minimap.
  --reference REFERENCE
                        fasta file with reference.
  --reads2contigs READS2CONTIGS
                        sam/bam file with reads mapped to contigs.
  -o OUTPUT, --output OUTPUT
                        output file (defaults to merged.sam)
  --mapped2ref MAPPED2REF
                        Read is mapped to reference. What should I do with it?
                        1. leave - leave it alone, don't touch it. 2. transfer
                        - if it can be transferred via contigs, do it.
  -k, --keep-unmapped-contigs
                        Add unmapped contigs to the references in the output?
                        By default I won't.
```

#### Mandatory arguments:

`-r` or `--reads` is the file with reads.
 If they are already mapped to the reference,
 it should be in bam or sam format.
 If they are not and you want IMCA to map them using mappy,
 they should be in fastq or fasta format.

`-c` or `--contigs` is the file with contigs.
 If they are already mapped to the reference,
 it should be in bam or sam format.
 If they are not and you want IMCA to map them using mappy,
 they should be in fastq or fasta format.
 To ensure better estimation of transfer coordinates,
 every mapped contig should have CIGAR string set.
 It is not compulsory, but it is highly encouraged.
 Many mappers will generate them by default
 (Bowtie2, minimap2 (if you request output in SAM format), BWA...).
 If you use some tool designed for aligning contigs to the reference (MUMmer, BLAST...)
 it might need some additional parsing to retrieve CIGAR strings
 from it's output.
 You can read more about how CIGAR string is constructed [here](https://www.drive5.com/usearch/manual/cigar.html).

##### Almost mandatory arguments:

Depending on what you provided as reads and contigs,
 you should also provide one of the arguments below.

`--reads2contigs` is the file with contigs mapped to the reference.
 It should be in bam or sam format.
 Provide it if you provided reads and contigs as sam/bam.
 To ensure better estimation of transfer coordinates,
 every mapped read should have CIGAR string set.
 It is not compulsory, but it is highly encouraged.

`--reference` is the sequence of reference genome.
 It should be in fasta or fastq format.
 Provide it if you provided reads and contigs as fasta/fastq.

#### Optional arguments:

`-o` is the desired named of the output file. Defaults to `merged.sam`.
 Currently output is always in .sam format.

`--mapped2ref` decides what to do with reads that are mapped both to the contigs and to the reference.
 Currently two options are available: `leave` means
 "leave them alone and keep the coordinates from `reads2reference` file"
 (that's the default behaviour),
 and `transfer` means "transfer coordinates whenever possible".

`-k` or `--keep-unmapped-contigs` is a flag that decides
 what to do with reads that are mapped to contigs which didn't map to the reference:
 should we just ignore and remove them
 (that's the default behaviour)
 or should we add these contigs as references to the output .sam file
 and keep the reads?
 The latter happens if you run the script with this option.
 Keep in mind that the header of your output bam file
 will now be different than the header of `reads2reference`;
 it will contain some additional references.

### Example datasets

To add.

### References

Heng Li. Minimap2: pairwise alignment for nucleotide sequences;
Bioinformatics, Volume 34, Issue 18, 15 September 2018, Pages 3094–3100, https://doi.org/10.1093/bioinformatics/bty191

D.R. Zerbino and E. Birney. 2008. Velvet: algorithms for de novo
short read assembly using de Bruijn graphs.
Genome Research,18:821-829


