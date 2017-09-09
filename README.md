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

* Assume you have downloaded the annotation and assembly of Rhodosporidium toruloides IFO0880_4 from JGI Mycocosm:
* Assume you have initialized a Pathway/Genome database named rhtocyc using Pathologic. `Tools->Pathologic->Database->Create new`

```
~/annot2pathologic/annot2pathologic.py \
   Mycocosm/Annotation/Filtered\ Models\ \(\"best\"\)/Genes/Rhoto_IFO0880_4_GeneCatalog_20170509.gff3 \
   Mycocosm/Annotation/Filtered\ Models\ \(\"best\"\)/Functional\ Annotations/EC\ annotations\ and\ KEGG\ pathway\ mappings/*.tab \
   Mycocosm/Annotation/Filtered\ Models\ \(\"best\"\)/Functional\ Annotations/KOG/Rhoto_IFO0880_4_GeneCatalog_proteins_20170509_KOG.tab \
   Mycocosm/Annotation/Filtered\ Models\ \(\"best\"\)/Functional\ Annotations/GO/Rhoto_IFO0880_4_GeneCatalog_proteins_20170509_GO.tab \
   Mycocosm/Assembly/Assembled_Scaffolds_Unmasked/Rhoto_IFO0880_4_AssemblyScaffolds.fasta \
   ~/ptools-local/pgdbs/user/rhtocyc/1.0/input
```


* You are now ready to go to `Build->Specify Replicons`


## Future work
* Generate an organism.dat file from the NCBI taxonomy that will allow you to run pathologic from the command-line, rather than invoking pathway-tools.
* Require only pointing to the base directory of a Mycoco
sm download, and annot2pathologic.py will find all the appropriate annotation files automagically.
* Require only pointing to the base URL of a Mycocosm download, and annot2pathologic.py will download the files it needs and convert them to pathologic file format automagically.