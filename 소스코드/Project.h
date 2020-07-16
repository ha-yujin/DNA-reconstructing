#pragma once
#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <time.h>
#include <fstream>
#include <map>
#define L 30 // 32~100
#define M 20000 // 2000�� ~ 2��
#define N 30000
#define D 2
using namespace std;

class Project
{
	char dna[4] = { 'A','C','G','T' };
	char *MyGenomeSeq = new char[N + 1];
	char *HumanSeq = new char[N + 1];
	string T$;
	string BWT_T = "";
	string *BWT_ta = new string[N + 1];
	int *sf = new int[N + 1]; // suffix array ����
	int **tally = new int*[N+1]; // tally table ����
	int last[4] = { 0, }; // last table ����
public:

	Project();
	~Project();
	void makeMyGenome();
	void makeShortRead();
	void makeHumanGenome();

	int charToInt(char a);

	// BWT algorithm //
	void myBWT();
	void makeTallyTab();
	void makeLastTab();
	void findPattern(int len);

	// �� algorithm //
	void trivial();
};

