
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <malloc.h>
#include <getopt.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>

#include"scaffoldSet.h"
#include "contigSet.h"
#include "aligningFromBam.h"
#include "fillGap.h"
using namespace std;


int main(int argc, char** argv) {

	if (argc == 1) {
		cout << "Please input correct command line.\nHRGF -s [scaffolds.fa] -r [hifi-reads.fa] -p [output-directory] -t [thread-count]\n";
		exit(0);
	}
	
	char* scaffoldSetFile = NULL;
	char* resultOutPutDirectory = NULL;
	char* longReadFile = NULL;
	char* shortReadFile1 = NULL;
	char* shortReadFile2 = NULL;
	long int threadCount = 1;


	int ch = 0;
	while ((ch = getopt(argc, argv, "s:r:p:t:h")) != -1) {
		switch (ch) {
		case 's': scaffoldSetFile = (char*)(optarg); break;
		case 'lr': longReadFile = (char*)(optarg); break;
		case 'sr1': shortReadFile1 = (char*)(optarg); break;
		case 'sr2': shortReadFile2 = (char*)(optarg); break;
		case 'p': resultOutPutDirectory = (char*)optarg; break;
		case 't': threadCount = atoi(optarg); break;
		case 'h': cout << "LSRGF -s [scaffolds.fasta] -lr [longReads.fasta] -sr1 [shortRead1.fasta] -sr2 [shortRead2.fasta] -p [output-directory] -t [thread-count]"; exit(0); break;
		case '?': cout << "Please input correct command line.\nHRGF -s [scaffolds.fa] -r [hifi-reads.fa] -p [output-directory] -t [thread-count]\n"; exit(0); break;
		case ':': cout << "Please input correct command line.\nHRGF -s [scaffolds.fa] -r [hifi-reads.fa] -p [output-directory] -t [thread-count]\n"; exit(0); break;
		default: cout << "Please input correct command line.\nHRGF -s [scaffolds.fa] -r [hifi-reads.fa] -p [output-directory] -t [thread-count]\n"; exit(0); break;

		}
	}

	if (opendir(resultOutPutDirectory) == NULL) {
		mkdir(resultOutPutDirectory, 0777);
	}

	char* longScaffoldSetFile = (char*)malloc(sizeof(char) * 150);
	strcpy(longScaffoldSetFile, resultOutPutDirectory);
	strcat(longScaffoldSetFile, "/longScaffold.fa");

	char* contigSetFile = (char*)malloc(sizeof(char) * 150);
	strcpy(contigSetFile, resultOutPutDirectory);
	strcat(contigSetFile, "/contigSet.fa");

	char* longAligningSam = (char*)malloc(sizeof(char) * 150);
	strcpy(longAligningSam, resultOutPutDirectory);
	strcat(longAligningSam, "/longAlingningResult.sam");
	char* longAligningBam = (char*)malloc(sizeof(char) * 150);
	strcpy(longAligningBam, resultOutPutDirectory);
	strcat(longAligningBam, "/longAligningResult.bam");
	
	char* shortAligningSam1 = (char*)malloc(sizeof(char) * 150);
	strcpy(shortAligningSam1, resultOutPutDirectory);
	strcat(shortAligningSam1, "/shortAlingningResult.sam");
	char* shortAligningBam1 = (char*)malloc(sizeof(char) * 150);
	strcpy(shortAligningBam1, resultOutPutDirectory);
	strcat(shortAligningBam1, "/shortAligningResult1.bam");
	
	char* shortAligningSam2 = (char*)malloc(sizeof(char) * 150);
	strcpy(shortAligningSam2, resultOutPutDirectory);
	strcat(shortAligningSam2, "/shortAlingningResult2.sam");
	char* shortAligningBam2 = (char*)malloc(sizeof(char) * 150);
	strcpy(shortAligningBam2, resultOutPutDirectory);
	strcat(shortAligningBam2, "/shortAligningResult2.bam");

	char* file = (char*)malloc(sizeof(char) * 150);
	strcpy(file, resultOutPutDirectory);
	strcat(file, "/contigPathInlongRead.fa");

	char* file1 = (char*)malloc(sizeof(char) * 150);
	strcpy(file1, resultOutPutDirectory);
	strcat(file1, "/optimizeContigPathInLongRead.fa");

	char* file2 = (char*)malloc(sizeof(char) * 150);
	strcpy(file2, resultOutPutDirectory);
	strcat(file2, "/supplementContigPathInLongRead.fa");

	FILE* fp;
	if ((fp = fopen(file1, "w")) == NULL) {
		printf("%s, does not exist!", file1);
		exit(0);
	}
	FILE* fp2;
	if ((fp2 = fopen(file2, "w")) == NULL) {
		printf("%s, does not exist!", file2);
		exit(0);
	}
	
	char* file4 = NULL;
	if (file4 == NULL) {
		file4 = (char*)malloc(sizeof(char) * fileNameLen + 50);
		strcpy(file4, resultOutPutDirectory);
		strcat(file4, "/originalPathInShortContig.fa");
	}

	char* file5 = (char*)malloc(sizeof(char) * fileNameLen + 50);
	strcpy(file5, resultOutPutDirectory);
	strcat(file5, "/optimaizePathInLongRead.fa");
	

	FILE* fp5;
	if ((fp5 = fopen(file5, "w")) == NULL) {
		printf("%s, does not exist!", file5);
		exit(0);
	}

	char* file6 = (char*)malloc(sizeof(char) * fileNameLen + 50);
	strcpy(file6, resultOutPutDirectory);
	strcat(file6, "/optimaizePathInShortContig.fa");


	FILE* fp6;
	if ((fp6 = fopen(file6, "w")) == NULL) {
		printf("%s, does not exist!", file6);
		exit(0);
	}
	
	DelectShortContig(scaffoldSetFile, longScaffoldSetFile);//去除长度小于1000的cotnig

	//ScaffoldSetHead* scaffoldSetHead = GetScaffoldSetFromScaffoldFile(scaffoldSetFile);//获取初始scaffold以及GAP位置信息
	ScaffoldSetHead* scaffoldSetHead = GetScaffoldSetFromScaffoldFile(longScaffoldSetFile);//获取去除长度小于1000的scaffold以及GAP位置信息

	GetContigSetFromScaffoldSetHead(scaffoldSetHead);//获取contig的位置信息

	OutPutContigSetInScaffoldSetHead(scaffoldSetHead, contigSetFile);

	int maxSize = 2000;
	char* line = (char*)malloc(sizeof(char) * maxSize);

	char command[3000];
	
	sprintf(command, "bwa index %s", contigSetFile,);
	system(command);
	sprintf(command, "bwa mem -a %s %s > %s", contigSetFile, shortReadFile1, shortAligningSam1);
	system(command);
	sprintf(command, "samtools view -Sb %s >%s", shortAligningSam1, shortAligningBam1);
	system(command);
	
	sprintf(command, "bwa mem -a %s %s > %s", contigSetFile, shortReadFile2, shortAligningSam2);
	system(command);
	sprintf(command, "samtools view -Sb %s >%s", shortAligningSam2, shortAligningBam2);
	system(command);
	
	sprintf(command, "bwa mem -t8 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y %s %s > %s", contigSetFile,longReadFile, longAligningSam);
	system(command);
	sprintf(command, "samtools view -Sb %s >%s", longAligningSam, longAligningBam);
	system(command);

	long int scaffoldIndex = -1;
	long int gapIndex = -1;

	ContigSetHead* contigSetHead = GetContigSetFromContigSetFile(contigSetFile);
	
	ReadMapPosition* readMapPosition = GetReadInformation(contigSetHead, shortAligningBam1, shortAligningBam2, longAligningBam, file1, file2, file4, insertSize, readType);
	
	GetAligningResult(scaffoldSetHead,contigSetHead, shortAligningBam1,shortAligningBam2,longAligningBam, file);
	
	SimpleResultHead* simpleResultHead = OptimaizeLongAlign(contigSetHead, file5);

	GetEstimatedGapLength(scaffoldSetHead, contigSetHead, file, line, maxSize);
	
	GetReadPoinstionInGap(contigSetHead, scaffoldSetHead, file, line, maxSize, fp);

	GetReadInGap(contigSetHead, scaffoldSetHead, longReadFile, file1, line, maxSize, resultOutPutDirectory);
	
	FillingGap(scaffoldSetHead, resultOutPutDirectory);

}