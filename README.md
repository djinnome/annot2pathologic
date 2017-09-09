# annot2pathologic

## Prerequisites

For this to work, you must install the following packages:

```
conda install -c bioconda biopython
conda install -c bioconda gffutils
conda install pandas numpy
```

## Installation

```
git clone https://github.com/djinnome/annot2pathologic.git
```

## Usage


```
annot2pathologic/annot2pathologic.py --help
usage: annot2pathologic.py [-h] gff ec kog go seq outputdir

Convert annotations to Pathologic file format

positional arguments:
  gff         input gff file name
  ec          input EC number annotations file name
  kog         input kog annotations file name
  go          input go annotations file name
  seq         input (unmasked) sequence file
  outputdir   output directory

optional arguments:
  -h, --help  show this help message and exit
```

## Example invocation

* Assume you have downloaded the annotation and assembly of Rhizopus delemar  99-880
* Assume you have initialized a Pathway/Genome database named rhizopuscyc using Pathologic. `Tools->Pathologic->Database->Create new`
* Download the selected files from [JGI Mycocosm](http://genome.jgi.doe.gov/pages/dynamicOrganismDownload.jsf?organism=Rhior3)

![JGI Mycocosm Rhior3](Rhior3.png "JGI Mycocosm Rhior3 Download")


```
~/annot2pathologic/annot2pathologic.py \
    Rhior3.filtered_proteins.BroadGene.gff3 \
    Rhior3_proteins_KEGG.tab Rhior3_proteins_KOG.tab \
    Rhior3_proteins_GO.tab Rhior3_scaffolds.fasta \
    ~/ptools-local/pgdbs/user/rhopuscyc/1.0/input
```


* You are now ready to go to `Build->Specify Replicons`


## Future work
* Generate an organism.dat file from the NCBI taxonomy that will allow you to run pathologic from the command-line, rather than invoking pathway-tools.
* Require only pointing to the base directory of a Mycocosm download, and annot2pathologic.py will find all the appropriate annotation files automagically.
* Require only pointing to the base URL of a Mycocosm download, and annot2pathologic.py will download the files it needs and convert them to pathologic file format automagically.