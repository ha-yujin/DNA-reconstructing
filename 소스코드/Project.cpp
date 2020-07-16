#include "Project.h"
#include <algorithm>
#define _CRT_SECURE_NO_WARNINGS
Project::Project() { // data 생성
	makeMyGenome();
	makeShortRead();
	makeHumanGenome();
}
Project::~Project() {}
void Project::makeMyGenome() {

	srand((unsigned)time(NULL)); // seed값으로 현재 시간을 입력받음

	for (int i = 0; i < N; i++) { // 랜덤으로 MyGenomeSeq 생성
		MyGenomeSeq[i] = dna[rand() % 4];
	}
	MyGenomeSeq[N] = 0;
	ofstream fout("input.txt"); // 파일에 저장
	if (fout.is_open()) {
		for (int i = 0; i < N; i++)
			fout << MyGenomeSeq[i];
		fout.close();
	}
}
void Project::makeShortRead() { // myGenomeSeq로부터 short_read 생성
	srand((unsigned)time(NULL));
	ofstream fout("reads.txt");
	int *randIndex = new int[M];
	int index;
	char **short_read = new char*[M]; // 길이가 L인 short read M개 동적할당
	for (int i = 0; i < M; ++i) {
		short_read[i] = new char[L];
	}
	for (int i = 0; i < M; i++) {
		index = rand() % (N - L+1);
		for (int k = 0; k < L; k++) // index 위치(랜덤으로 지정)에서부터 L개의 염기서열 자름
			short_read[i][k] = MyGenomeSeq[index++];
	}
	if (fout.is_open()) { // 파일이 정상적으로 열린 경우
		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < L; j++)
				fout << short_read[i][j]; // 오픈한 파일에 내용 출력
			fout << endl;
		}
		fout.close(); // 파일을 닫아준다
	}
	for (int i = 0; i < M; i++) // 동적할당 해제
		delete[] short_read[i];
	delete[] short_read;
}
void Project::makeHumanGenome() { // myGenome으로부터 reference인 HumanGenome 생성
	srand((unsigned)time(NULL));
	int num = N*0.01; // 스닙 개수 MyGenome길이의 1%
	for (int i = 0; i < N; i++)
		HumanSeq[i] = MyGenomeSeq[i];
	for (int i = 0; i < num; i++) // 랜덤한 위치에 랜덤한 염기로 바꿔줌
	{
		int index = rand() % N;
		HumanSeq[index] = dna[rand() % 4];
	}
	HumanSeq[N] = 0;
	ofstream fout("humanSeq.txt"); // 파일에 출력
	if (fout.is_open()) {
		for (int i = 0; i < N; i++)
			fout << HumanSeq[i];
		fout.close();
	}
}
int Project::charToInt(char a) { // 문자를 숫자로 바꾼다
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
	
	/* BWT table 생성 */
	int i, j;
	BWT_ta[0] = T$[len - 1];
	for (i = len - 2; i >= 0; i--) { // ACGT일 때 --> T , GT , CGT, ACGT
		j = i;
		BWT_ta[len-1-i] = T$[j] + BWT_ta[len-i - 2];
	}
    sort(BWT_ta, BWT_ta+len); // BWT table 정렬
	for (int i = 0; i < len; i++) // suffix array 생성
		sf[i] = len - BWT_ta[i].size();

	/*for (int i= 0; i < len; i++) {
		cout << sf[i] << "  ";
		cout << BWT_ta[i] << endl;
	}*/
	for (int i = 0; i <len; i++) // BWT(T) 구하기
	{
		int idx = sf[i];
		if (idx == 0)
			BWT_T+= T$[len - 1];
		else
			BWT_T+= T$[sf[i]-1];
	}
	makeTallyTab(); // tally table 만들기
	makeLastTab(); 	//last table 만들기

	findPattern(len); // 만들어둔 정보들을 갖고 originalSeq 찾기
	for (int i = 0; i < N + 1; i++) // data들 동적할당 해제
		delete tally[i];
	delete[] tally;
	delete[] sf;
	delete[] BWT_ta;
}

/* 
 * last 배열과 함께 사용해 
 * L과 F를 매칭해주는 것을 수월하게 하기 위한 tally table 생성   
 * 문자의 존재 여부도 알 수 있음
 */
void Project::makeTallyTab() { 
	for (int i = 0; i < N + 1; i++) // 동적 할당
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
		if (i != 0) // tally table은 누적 table이므로 전 행과 더해줌
		{
			for (int j = 0; j < 4; j++)
				tally[i][j] += tally[i - 1][j];
		}
	}
}
/*
 * 이전에 생성한 tally table을 갖고 last table을 생성하는 함수
 * 패턴을 찾을 때 suffix array에서 특정 문자의 시작과 끝을 찾기 위함
 */
void Project::makeLastTab() {
	last[0] = tally[N][0]; // 맨 마지막 행의 tally table을 이용하면 last table 쉽게 생성 가능
	for (int i = 1; i <4; i++)
		last[i] = tally[N][i] + last[i - 1];
}
/* read와 humanSeq를 비교하면서 originalSeq를 찾아나가는 과정 */
void Project::findPattern(int len) {

	char *findOriginal = new char[N + 1]();
	ifstream fin("reads.txt");
	string reads;
	map<int, int> tmp; // 미스매치 저장위한 변수
	int start, end, cur; // 문자의 시작과 끝, 현재 위치 가리키는 변수
	int i, cnt = 0, mismatch=0;
	int ch, nextch;
	while (fin.peek() != EOF) { // 파일의 끝까지 읽어온다
		getline(fin, reads);
		int refIn = L - 1; 
		ch = charToInt(reads[refIn--]); // read의 맨 마지막 글자를 숫자로 변환

		/* 각 문자에 따라 F에서 시작, 끝 범위 지정*/
		if (ch == 0)
			start = 1;
		else
			start = last[ch - 1] + 1;
		end = last[ch];
		for (i = start; i < end; i++) { 
			cur = i;
			refIn = L - 2;
			while (refIn>=0) { // 한 read를 다 읽을 때까지 반복
				if (BWT_T[cur] != reads[refIn]) // read 문자와 일치하지 않을 때
				{
					mismatch++; // mismatch 증가
					if (mismatch > D) // 정해준 수치를 넘어갈 경우 탐색 종료
					{
						mismatch = 0;
						break;
					}
				}
				refIn--;

				/* 다음 문자 위치를 찾는 과정 , 만약 $만날 경우 탐색 종료 */
				nextch = charToInt(BWT_T[cur]);
				if (nextch == -1)
					break;
				else {
					if (nextch == 0) // A인 경우
						cur = tally[cur][nextch];
					else
						cur = tally[cur][nextch] + last[nextch - 1]; // C,G,T인 경우
				}
			}
			if (refIn < 0) // 한 read 탐색이 완료된 경우
			{
				tmp[mismatch] = sf[cur]; // mismatch를 key로 index를 value로 저장
				if (mismatch == 0) // 만약 정확한 위치를 찾은 경우 탐색 종료
					break;
			}
			mismatch = 0;
		}
		int idx;

		/* mismatch가 가장 적은 인덱스를 read의 인덱스로 보고 findOriginal에 넣어준다 */
		if (tmp.count(0) == 0)
		{
			if (tmp.count(1) == 0)
			{
				if (tmp.count(2) == 0)
					idx = tmp[3]; // mismatch 3개 허용할 경우
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
		tmp.clear(); mismatch = 0; // mismatch 관련 변수 초기화
	}

	for (int j = 0; j < N; j++) { 
		if (findOriginal[j] == 0) // cover가 되지 않은 부분은 humanSeq로 채운다
			findOriginal[j] = HumanSeq[j];
		if (findOriginal[j] == MyGenomeSeq[j]) // 정확도 계산
			cnt++;
	}findOriginal[N] = 0;
	cout <<"BWT 정확도 : " << (double)cnt / (double)N * 100 <<"%"<< endl;
	delete[] findOriginal;
}
void Project::trivial() {

	ifstream fin("reads.txt");
	string reads;

	
	char *findOriginal = new char[N + 1]();
	int i, j;
	int cnt = 0;
	int mismatch = 0;
	map<int, int> tmp; // 미스매치 저장위한 변수
	if (fin.is_open())
	{
		while (fin.peek() != EOF) { // peek()로 스트림에 뭐가 있는지 살표봄
			getline(fin, reads); // 입력스트림에서 한줄 읽어옴

			for (i = 0; i <= N - L; i++) { // 길이가 N인 humanSeq 탐색
				for (j = 0; j < L; j++) // 한 read를 탐색
				{
					if (HumanSeq[i + j] != reads[j]) // 문자가 같지 않을 경우 mismatch 증가
					{
						mismatch++;
					}
				}
				if (j == L) // 한 read 탐색이 끝난 경우
				{
					if (mismatch <= D)
						tmp[mismatch] = i;
					mismatch = 0;
				}
			}

			/* mismatch 개수가 가장 적은 인덱스를 read의 위치로 본다 */
			if (tmp.count(0) == 0)
			{
				if (tmp.count(1) == 0)
				{
					if (tmp.count(2) == 0)
						i = tmp[3]; // mismatch 3개 허용할 경우
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
			tmp.clear(); // 변수 초기화
		}findOriginal[N] = 0;


		for (int i = 0; i < N; i++) // 추측해서 만든 findOriginal과 MyGenomeSeq비교
		{
			if (findOriginal[i] == 0)
				findOriginal[i] = HumanSeq[i];
			if (findOriginal[i] == MyGenomeSeq[i])
				cnt++;
		}
		cout << "Trivial 정확도 : " << (double)cnt / (double)N * 100 <<"%"<< endl;

		delete[] findOriginal;
		fin.close();
	}
}