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
	cout << "BWT����ð� : "<<(double)(end - begin) / (double)CLOCKS_PER_SEC <<"��"<< endl;
	
	cout << endl << "---------------------" << endl;
	begin = clock();
	pro.trivial();
	end = clock();
	cout << "Trivial����ð� : "<<(double)(end - begin) / (double)CLOCKS_PER_SEC <<"��"<< endl;

}