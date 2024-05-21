# LSRGF
LSRGF is a Hybrid gap filling based on short and long reads
=========


gap-filling: LSEGF
=================

1) Introduction
```
    LSRGF is an gap-filling tool which aims to filling gaps of scaffolds. 
    The input data of LSRGF is the short reads (fasta format),long reads (fasta format) and the scaffolds (fasta format). 
```
2) Before installing and running
```
    Please install BWA from https://github.com/lh3/bwa.
	Please install Samtools from https://sourceforge.net/projects/samtools/files/samtools/.
	Please build and install Bamtools from https://github.com/pezmaster31/bamtools.
```
3) Installing.
```
    LSRGF should run on Linux operating sysetm with gcc. We test LSRGF using gcc9.4.0 on Ubuntu.
    Create a main directory. Copy all source code to this directory.
	cd LSRGF 
	export BAMTOOLS_HOME_INCLUDE=/path_bamtools_include_api_shared/
	export BAMTOOLS_HOME_LIB=/path_bamtools_lib_libbamtools.a/
	make all
```
4) Running.
```
	LSRGF -s [scaffolds.fasta] -lr [longReads.fasta] -sr1 [shortRead1.fasta] -sr2 [shortRead2.fasta] -p [output-directory] -t [thread-count]
    
    -s <scaffolds.fasta>: 
	    The file includes scaffolds produced by one assembler.
	-lr <longReadType.fasta>: 
	    The type of long reads.
	-sr1 <shortReadeadType1.fasta>: 
	    The type of short reads.
	-sr2 <shortReadType2.fasta>: 
	    The type of short reads.	
	-p <outdir>: 
	   The output path of the file.
	-t <threadNumber>: 
	   The number of threads used by the tool.
```
5) Output.
```
    The output file "gapfilling-scaffoldset.fa" is the gap-filling result. 
```


