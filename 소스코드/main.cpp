#include "Project.h"
int main() {

	Project pro;
	clock_t begin, end;
	cout << "L=" << L << endl;
	cout << "M=" << M<<endl;
	cout << "N=" << N<<endl;
	cout << "---------------------" << endl;

	begin = clock();
	pro.myBWT();
	end = clock();
	cout << "BWT실행시간 : "<<(double)(end - begin) / (double)CLOCKS_PER_SEC <<"초"<< endl;
	
	cout << endl << "---------------------" << endl;
	begin = clock();
	pro.trivial();
	end = clock();
	cout << "Trivial실행시간 : "<<(double)(end - begin) / (double)CLOCKS_PER_SEC <<"초"<< endl;

}