# How to perform an initial QC of your sequencing data?

Enter the data directory and make yourself familiar with the files.
Usually there should be one or two directories containing all your .fastq/.fq sequencing files (more information about this format [here](https://en.wikipedia.org/wiki/FASTQ_format)). In most cases the files are compressed (.gz or .bz). Additionally, there should be some kind of Readme file with specific information about the received data and eventually a design file, translating file names into your sample names. In some cases the sequencing agency already performed some cleaning procedure. In this case you most likely received a directory with the data before and after cleaning, respectively. Raw (and clean) data can be treated by different filtering and trimming steps:

## Check, if you downloaded the files correctly

Some sequencing agencies send md5 checksums to check data integrity after downloading your data, e.g. the md5 checksum should be identical before and after the download. Otherwise your downloaded file might be corrupt.

Look for a file called 'md5.txt' or similar and run the following command to generate an output file called md5.check:

```bash
#change "md5.txt" to the appropiate file name
md5sum -c md5.txt > md5.check
```

If your report md5.check shows something like this you can continue with your QC analysis:

```bash
some_name_16S_205_1.fq.gz: OK
some_name_16S_205_2.fq.gz: OK
some_othername_16S_205_1.fq.gz: OK
some_othername_16S_205_2.fq.gz: OK
```

If your report shows something like this you should repeat the download of the specific file(s) and if the error still occurs, you should contact your sequencing facility:

```bash
some_name_16S_205_1.fq.gz: FAILED
some_name_16S_205_2.fq.gz: OK
some_othername_16S_205_1.fq.gz: OK
some_othername_16S_205_2.fq.gz: OK
```

## Run FastQC and Fast-screen to generate various statistics about your sequencing reads

Clone the gitlab repository with the required scripts and directories. If you are a new user, then first log in with your NIOO credentials on Gitlab (https://gitlab.bioinf.nioo.knaw.nl). This will activate your gitlab account and will allow you to clone repositories.

```bash
# replace <user-name> by your own NIOO login-name, e.g FleurG
git clone https://gitlab.bioinf.nioo.knaw.nl/tutorials/QC-NGS.git
```

Enter the QC-NGS directory and create a conda environment containing all the required tools and dependencies.

```bash
conda env create -f environment.yml
# to inactivate the conda environment type:
# (e.g. when you finished all your analysis)
source deactivate
# to connect to your conda environment again type:
# (e.g. the next time you want to perform your QC analysis)
source activate QC
# if you do not end up in (QC) but get an error
# or you get into (qiime), then first execute:
export PATH=/opt/miniconda3/bin:$PATH
```

##### Analyse your raw (and cleaned) sequencing reads with FastQC.

You can run several files in one directory at once by using the asterix symbol. By default the output will be stored in the sub-directories `./reports/clean` and `./reports/raw.` We here expect that you have a data directory with raw and a directory with cleaned reads and that you linked them to the directories ./data/clean and ./data/raw, respectively. Use `ln -s`.

```bash
# please adjust the path to your reads, accordingly
fastqc -o ./reports/clean/ -t 12 ./data/clean/*.fastq.gz
# please ajdust the path to your reads, accordingly
fastqc -o ./reports/raw/ -t 12 ./data/raw/*.fastq.gz
```

Move on with the tutorial by executing Fastq-screen. You can scip the section about making a fastq-screen database, when using the pre-existing databases is sufficient.

#####  Make a fastq-screen database

Fast-screen will use bowtie2 by default to align your data. Optional you can use BWA or bowtie by adding the --aligner flag.

At the moment Fast-screen will screen the following sources of possible contaminations:
- vectors
- adapters
- phix
- Escherichia coli
- Homo sapiens
- Arabidopsis thaliana

If you miss anything in this list, please contact BU bioinformaticians and add the desired database by:

- Downloading a fasta-file with your sequence(s)and indexing your database/genome. By default, please index your genome/database with bowtie2.

```bash
# For bowtie index the database with:
bowtie-build sequence.fasta sequence
bowtie-inspect sequence
# For bowtie2 index the database with:
bowtie2-build sequence.fasta sequence
bowtie2-inspect sequence
# For bwa index the database
bwa index -p sequence ./sequence.fasta
```

- Adding the following text at the bottom of the fastq_screen.conf file:

```bash
## Sequence - Some information about the source of the sequence
## A link to the url, where you downloaded the sequence
DATABASE        Sequence     /path/to/your/indexed/sequence
```

##### Start fastq-screen to identify contaminants in your raw and (cleaned sequencing) reads

```bash
cd ./fastq_screen_v0.11.1
./fastq_screen --outdir ../reports/clean-screen/ ../data/clean/*.fastq.gz #add --aligner bowtie or --aligner bwa, if you do not want to use bowtie (default)
./fastq_screen --outdir ../reports/raw-screen/ ../data/raw/*.fastq.gz #add --aligner bowtie or --aligner bwa, if you do not want to use bowtie (default)
```

##### Summarize the reports with multiqc

```bash
cd ..
multiqc -d -n report_all ./reports
```

## How to read and interpret the final report?
- Download the report_all.html file from the server and open it in your web-browser.
- You can watch the tutorial video as indicated in the report. The multiqc report allows you to rename, filter or select samples and export the plots or the data.
- Multiqc collects 'general statistics' about the number of duplicated reads (%Dups), the GC content (%GC) and the number of sequencing reads (M Seqs). This can give you a quick overlook about your data.
- The output of fastq-screen is visualized in a bar chart. For each sample the percentage of reads that matched with a "contaminants" database is plotted. Of course you have to put the results into the right context: if you sequenced 16S amplicons from bacteria, you will most likely get a high percentage of hits with the E. coli genome.
- The FastQC report consists of 9 different parts and will indicate the failed and passed samples. You can access a detailed description and explanation about each part in the [FastQC online help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help) Also here it is important to evaluate the results in the biological context, e.g. with amplicon sequencing you expect a bias in GC content and a high duplication rate:
  - Sequence quality histograms: This shows an overview of the range of quality values across all bases at each position in the FastQ file. It is normal that the quality values decrease at higher base positions. You can think about trimming your data based on quality if an overall high base quality is required for your analysis.
  - Per sequence quality scores: The number of reads is plotted over the average quality per read. It is often the case that a subset of sequences will have universally poor quality, often because they are poorly imaged (on the edge of the field of view etc), however these should represent only a small percentage of the total sequences.
  - Per base sequence content: This plots out the proportion of each base position in a file for which each of the four normal DNA bases has been called. In a random library you would expect that there would be little to no difference between the different bases of a sequence run, so the lines in this plot should run parallel with each other. The relative amount of each base should reflect the overall amount of these bases in your genome, but in any case they should not be hugely imbalanced from each other.
  - Per sequence GC content: This measures the GC content across the whole length of each sequence in a file and compares it to a modeled normal distribution of GC content. In a normal random library you would expect to see a roughly normal distribution of GC content where the central peak corresponds to the overall GC content of the underlying genome.
  - Per base N content: This plots out the percentage of base calls at each position for which an N was called. If you are unsatisfied with the result you can consider filtering reads that contain one (two, three...) or more N's.
  - Sequence length distribution: This shows the distribution of fragment sizes in the file which was analysed. If all reads have the same length, than the length of the reads is shown instead.
  - Duplicate sequences: This counts the degree of duplication for every sequence in a library and creates a plot showing the relative number of sequences with different degrees of duplication. If you encounter high duplication reads and did not execute a targeted study like amplicon sequencing, your data might be affected by a PCR bias.
  - Overrepresented sequences: This plots the total proportion of sequences which make up more than 0.1% of the total.
  - Adapter content: The plot shows a cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.

## How to proceed if quality should be improved?

##### Demultiplexing
Samples can be multiplexed to run in a single high-throughput sequencing instrument. To accomplish this, individual "barcode" sequences are added to each sample so they can be distinguished and sorted during data analysis. Pooling samples exponentially increases the number of samples analyzed in a single run, without drastically increasing cost or time. Usually, your sequencing agency will provide you with the demultiplexed data. If not, you can demultiplex the reads using [FASTX](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_barcode_splitter_usage), [IDEMP](https://github.com/yhwu/idemp) or any other tool that you prefer

##### Quality trimming
Low-quality base calls are trimmed off from the 3' end of the reads. This efficiently removes poor quality portions of the reads. You can use [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming), TrimGalore! or [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

##### Adapter trimming
To execute a sequencing run, adapters are ligated to your DNA-fragments. Most of the times, your sequencing agency removes those sequences by default. If not, you can performing adapter trimming with tools like [cutadapt](http://cutadapt.readthedocs.io/en/stable/index.html), [FASTX](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_clipper_usage) or [TrimGalore!](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)

##### Length trimming
For some applications it is required to trim all reads to the same length. This can be done using [FASTX](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_trimmer_usage).

##### Length filtering
Quality and/or adapter trimming may result in very short sequences (sometimes as short as 0 bp). You can filter trimmed reads based on their sequence length. This is to reduce the size of the output file and to avoid crashes of alignment programs which require sequences with a certain minimum length. You can use [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html#filtering-reads) or TrimGalore!.

##### Quality filtering
You can filter reads based on overall read-quality. You can use [FASTX](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastq_quality_filter_usage).

##### More programs

- seqtk
- fastp
- Trimmomatic
