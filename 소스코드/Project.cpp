#include "Project.h"
#include <algorithm>
#define _CRT_SECURE_NO_WARNINGS
Project::Project() { // data ����
	makeMyGenome();
	makeShortRead();
	makeHumanGenome();
}
Project::~Project() {}
void Project::makeMyGenome() {

	srand((unsigned)time(NULL)); // seed������ ���� �ð��� �Է¹���

	for (int i = 0; i < N; i++) { // �������� MyGenomeSeq ����
		MyGenomeSeq[i] = dna[rand() % 4];
	}
	MyGenomeSeq[N] = 0;
	ofstream fout("input.txt"); // ���Ͽ� ����
	if (fout.is_open()) {
		for (int i = 0; i < N; i++)
			fout << MyGenomeSeq[i];
		fout.close();
	}
}
void Project::makeShortRead() { // myGenomeSeq�κ��� short_read ����
	srand((unsigned)time(NULL));
	ofstream fout("reads.txt");
	int *randIndex = new int[M];
	int index;
	char **short_read = new char*[M]; // ���̰� L�� short read M�� �����Ҵ�
	for (int i = 0; i < M; ++i) {
		short_read[i] = new char[L];
	}
	for (int i = 0; i < M; i++) {
		index = rand() % (N - L+1);
		for (int k = 0; k < L; k++) // index ��ġ(�������� ����)�������� L���� ���⼭�� �ڸ�
			short_read[i][k] = MyGenomeSeq[index++];
	}
	if (fout.is_open()) { // ������ ���������� ���� ���
		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < L; j++)
				fout << short_read[i][j]; // ������ ���Ͽ� ���� ���
			fout << endl;
		}
		fout.close(); // ������ �ݾ��ش�
	}
	for (int i = 0; i < M; i++) // �����Ҵ� ����
		delete[] short_read[i];
	delete[] short_read;
}
void Project::makeHumanGenome() { // myGenome���κ��� reference�� HumanGenome ����
	srand((unsigned)time(NULL));
	int num = N*0.01; // ���� ���� MyGenome������ 1%
	for (int i = 0; i < N; i++)
		HumanSeq[i] = MyGenomeSeq[i];
	for (int i = 0; i < num; i++) // ������ ��ġ�� ������ ����� �ٲ���
	{
		int index = rand() % N;
		HumanSeq[index] = dna[rand() % 4];
	}
	HumanSeq[N] = 0;
	ofstream fout("humanSeq.txt"); // ���Ͽ� ���
	if (fout.is_open()) {
		for (int i = 0; i < N; i++)
			fout << HumanSeq[i];
		fout.close();
	}
}
int Project::charToInt(char a) { // ���ڸ� ���ڷ� �ٲ۴�
	if (a == 'A')
		return 0;
	else if (a == 'C')
		return 1;
	else if (a == 'G')
		return 2;
	else if (a == 'T')
		return 3;
	else
		return -1;
}
void Project::myBWT() {

	T$ = HumanSeq;
	T$ += "$";
	int len = N + 1;
	
	/* BWT table ���� */
	int i, j;
	BWT_ta[0] = T$[len - 1];
	for (i = len - 2; i >= 0; i--) { // ACGT�� �� --> T , GT , CGT, ACGT
		j = i;
		BWT_ta[len-1-i] = T$[j] + BWT_ta[len-i - 2];
	}
    sort(BWT_ta, BWT_ta+len); // BWT table ����
	for (int i = 0; i < len; i++) // suffix array ����
		sf[i] = len - BWT_ta[i].size();

	/*for (int i= 0; i < len; i++) {
		cout << sf[i] << "  ";
		cout << BWT_ta[i] << endl;
	}*/
	for (int i = 0; i <len; i++) // BWT(T) ���ϱ�
	{
		int idx = sf[i];
		if (idx == 0)
			BWT_T+= T$[len - 1];
		else
			BWT_T+= T$[sf[i]-1];
	}
	makeTallyTab(); // tally table �����
	makeLastTab(); 	//last table �����

	findPattern(len); // ������ �������� ���� originalSeq ã��
	for (int i = 0; i < N + 1; i++) // data�� �����Ҵ� ����
		delete tally[i];
	delete[] tally;
	delete[] sf;
	delete[] BWT_ta;
}

/* 
 * last �迭�� �Բ� ����� 
 * L�� F�� ��Ī���ִ� ���� �����ϰ� �ϱ� ���� tally table ����   
 * ������ ���� ���ε� �� �� ����
 */
void Project::makeTallyTab() { 
	for (int i = 0; i < N + 1; i++) // ���� �Ҵ�
		tally[i] = new int[4]();
	for (int i = 0; i < N+ 1; i++)
	{
		if (BWT_T[i] == 'A')
			tally[i][0]++;
		else if (BWT_T[i] == 'C')
			tally[i][1]++;
		else if (BWT_T[i] == 'G')
			tally[i][2]++;
		else if (BWT_T[i] == 'T')
			tally[i][3]++;
		if (i != 0) // tally table�� ���� table�̹Ƿ� �� ��� ������
		{
			for (int j = 0; j < 4; j++)
				tally[i][j] += tally[i - 1][j];
		}
	}
}
/*
 * ������ ������ tally table�� ���� last table�� �����ϴ� �Լ�
 * ������ ã�� �� suffix array���� Ư�� ������ ���۰� ���� ã�� ����
 */
void Project::makeLastTab() {
	last[0] = tally[N][0]; // �� ������ ���� tally table�� �̿��ϸ� last table ���� ���� ����
	for (int i = 1; i <4; i++)
		last[i] = tally[N][i] + last[i - 1];
}
/* read�� humanSeq�� ���ϸ鼭 originalSeq�� ã�Ƴ����� ���� */
void Project::findPattern(int len) {

	char *findOriginal = new char[N + 1]();
	ifstream fin("reads.txt");
	string reads;
	map<int, int> tmp; // �̽���ġ �������� ����
	int start, end, cur; // ������ ���۰� ��, ���� ��ġ ����Ű�� ����
	int i, cnt = 0, mismatch=0;
	int ch, nextch;
	while (fin.peek() != EOF) { // ������ ������ �о�´�
		getline(fin, reads);
		int refIn = L - 1; 
		ch = charToInt(reads[refIn--]); // read�� �� ������ ���ڸ� ���ڷ� ��ȯ

		/* �� ���ڿ� ���� F���� ����, �� ���� ����*/
		if (ch == 0)
			start = 1;
		else
			start = last[ch - 1] + 1;
		end = last[ch];
		for (i = start; i < end; i++) { 
			cur = i;
			refIn = L - 2;
			while (refIn>=0) { // �� read�� �� ���� ������ �ݺ�
				if (BWT_T[cur] != reads[refIn]) // read ���ڿ� ��ġ���� ���� ��
				{
					mismatch++; // mismatch ����
					if (mismatch > D) // ������ ��ġ�� �Ѿ ��� Ž�� ����
					{
						mismatch = 0;
						break;
					}
				}
				refIn--;

				/* ���� ���� ��ġ�� ã�� ���� , ���� $���� ��� Ž�� ���� */
				nextch = charToInt(BWT_T[cur]);
				if (nextch == -1)
					break;
				else {
					if (nextch == 0) // A�� ���
						cur = tally[cur][nextch];
					else
						cur = tally[cur][nextch] + last[nextch - 1]; // C,G,T�� ���
				}
			}
			if (refIn < 0) // �� read Ž���� �Ϸ�� ���
			{
				tmp[mismatch] = sf[cur]; // mismatch�� key�� index�� value�� ����
				if (mismatch == 0) // ���� ��Ȯ�� ��ġ�� ã�� ��� Ž�� ����
					break;
			}
			mismatch = 0;
		}
		int idx;

		/* mismatch�� ���� ���� �ε����� read�� �ε����� ���� findOriginal�� �־��ش� */
		if (tmp.count(0) == 0)
		{
			if (tmp.count(1) == 0)
			{
				if (tmp.count(2) == 0)
					idx = tmp[3]; // mismatch 3�� ����� ���
				else
					idx = tmp[2];
			}
			else
				idx = tmp[1];
		}
		else
			idx = tmp[0];
		for (int i = 0; i < L; i++) {
			findOriginal[idx + i] = reads[i];
		}
		tmp.clear(); mismatch = 0; // mismatch ���� ���� �ʱ�ȭ
	}

	for (int j = 0; j < N; j++) { 
		if (findOriginal[j] == 0) // cover�� ���� ���� �κ��� humanSeq�� ä���
			findOriginal[j] = HumanSeq[j];
		if (findOriginal[j] == MyGenomeSeq[j]) // ��Ȯ�� ���
			cnt++;
	}findOriginal[N] = 0;
	cout <<"BWT ��Ȯ�� : " << (double)cnt / (double)N * 100 <<"%"<< endl;
	delete[] findOriginal;
}
void Project::trivial() {

	ifstream fin("reads.txt");
	string reads;

	
	char *findOriginal = new char[N + 1]();
	int i, j;
	int cnt = 0;
	int mismatch = 0;
	map<int, int> tmp; // �̽���ġ �������� ����
	if (fin.is_open())
	{
		while (fin.peek() != EOF) { // peek()�� ��Ʈ���� ���� �ִ��� ��ǥ��
			getline(fin, reads); // �Է½�Ʈ������ ���� �о��

			for (i = 0; i <= N - L; i++) { // ���̰� N�� humanSeq Ž��
				for (j = 0; j < L; j++) // �� read�� Ž��
				{
					if (HumanSeq[i + j] != reads[j]) // ���ڰ� ���� ���� ��� mismatch ����
					{
						mismatch++;
					}
				}
				if (j == L) // �� read Ž���� ���� ���
				{
					if (mismatch <= D)
						tmp[mismatch] = i;
					mismatch = 0;
				}
			}

			/* mismatch ������ ���� ���� �ε����� read�� ��ġ�� ���� */
			if (tmp.count(0) == 0)
			{
				if (tmp.count(1) == 0)
				{
					if (tmp.count(2) == 0)
						i = tmp[3]; // mismatch 3�� ����� ���
					else
						i = tmp[2];
				}
				else
					i = tmp[1];
			}
			else
				i = tmp[0];

			for (int k = 0; k < L; k++) {
				findOriginal[i + k] = reads[k];
			}
			tmp.clear(); // ���� �ʱ�ȭ
		}findOriginal[N] = 0;


		for (int i = 0; i < N; i++) // �����ؼ� ���� findOriginal�� MyGenomeSeq��
		{
			if (findOriginal[i] == 0)
				findOriginal[i] = HumanSeq[i];
			if (findOriginal[i] == MyGenomeSeq[i])
				cnt++;
		}
		cout << "Trivial ��Ȯ�� : " << (double)cnt / (double)N * 100 <<"%"<< endl;

		delete[] findOriginal;
		fin.close();
	}
}