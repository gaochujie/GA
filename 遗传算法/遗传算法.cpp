// 遗传算法.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include<algorithm>
using namespace std;

const int N = 100;
double Pc = 0.7;
double Pm = 0.07;
const int MaxGens = 500;
const int Vars = 2;
const double PI = 3.1415926;
double avg[10];
int iii = 0;
struct genetype{     //染色体
	double gene[Vars];   //基因
	double fitness;    //适应值
	double pi;   //选择概率
	double qi;   //累积概率
	//double lval[Vars];   //左值
	//double rval[Vars];  //右值
};
class GA {
public:
	genetype population[N + 1];
	genetype newpopulation[N + 1];
	const double x[Vars][2] = { {-3.0,12.1},{4.1,5.8} };
	GA(){}
	void initialize();   //初始化
	void countfit();    //计算适应值
	void countmaxfit();   //计算最优解
	double fx(double x1,double x2);       //目标函数
	double randval(double low, double high);  //随机数生成函数
	void RWS();      //轮盘赌选择算法
	void crossover();   //交配
	void Xover(int,int);   //交换部分基因
	void swap(double*, double*);  //交换
	void variation();   //变异
	void printval();    //输出
};

double GA::randval(double low, double high) {
	double val;
	val = low + (high - low)*rand()*1.0 / RAND_MAX;
	return val;
}

void GA::initialize() {
	for (int i = 0;i < N;i++) {
		for (int j = 0;j < Vars;j++) {
			population[i].gene[j] = randval(x[j][0], x[j][1]);
		}
		double x1 = population[i].gene[0];
		double x2 = population[i].gene[1];
		population[i].fitness = fx(x1, x2);
	}
	population[N] = population[N - 1];
}

double GA::fx(double x1, double x2) {
	double y = (21.5 + x1 * sin(4 * PI*x1) + x2*sin(20 * PI*x2));
	return y;
}

void GA::countfit() {
	for (int i = 0;i < N;i++) {
		double x1 = population[i].gene[0];
		double x2 = population[i].gene[1];
		population[i].fitness = fx(x1, x2);
	}
}
void GA::countmaxfit() {
	genetype MAX;
	MAX.fitness = population[0].fitness;
	for (int i = 1;i < N;i++) {
		if (population[i].fitness >= MAX.fitness) {
			MAX.fitness = population[i].fitness;
			for (int j = 0;j < Vars;j++) {
				MAX.gene[j] = population[i].gene[j];
			}
		}
	}
	if (population[N].fitness < MAX.fitness) {
		population[N].fitness = MAX.fitness;
		for (int i = 0;i < Vars;i++)
			population[N].gene[i] = MAX.gene[i];
	}
}
void GA::RWS() {
	double P = 0.0;
	for (int i = 0;i < N;i++) {
		P += population[i].fitness;
	}
	population[0].pi = population[0].fitness / P;
	population[0].qi = population[0].pi;
	for (int i = 1;i < N;i++) {
		population[i].pi = population[i].fitness / P;
		population[i].qi = population[i].pi + population[i - 1].qi;
		//cout << population[i].qi << endl;
	}
	for (int i = 0;i < N;i++) {
		double r = randval(0.0, 1.0);
		for (int j = 0;j < N;j++) {
			if (r <= population[j].qi) {
				newpopulation[i] = population[j];
				break;
			}
		}
	}
	newpopulation[N] = population[N];
	for (int i = 0;i < N;i++) {
		population[i] = newpopulation[i];
	}
	//for (int i = 0;i < N+1;i++)
	//	cout << newpopulation[i].fitness<<endl;
	//cout << randval(0.0, 1.0) << endl;
}
void GA::swap(double *x, double *y) {
	double temp;
	temp = *x;
	*x = *y;
	*y = temp;
}
void GA::Xover(int one,int two) {
	int point = 0;
	if (Vars > 1) {
		if (Vars == 2) {
			point = 1;
		}
		else {
			point = (rand() % (Vars - 1) + 1);
		}
		for (int i = 0;i < point;i++) {
			swap(&population[one].gene[i], &population[two].gene[i]);
		}
	}
}
void GA::crossover() {
	int first = 0;
	int one = 0;
	double x = 0.0;
	for (int i = 0;i < N;i++) {
		x = randval(0.0, 1.0);
		if (x < Pc) {
			first++;
			if (first % 2 == 0) {
				Xover(one, i);
			}
			else {
				one = i;
			}
		}
	}
}

void GA::variation() {
	double r = 0.0;
	for (int i = 0;i < N;i++) {
		for (int j = 0;j < Vars;j++) {
			r = randval(0.0, 1.0);
			if (r < Pm) {
				population[i].gene[j] = randval(x[j][0], x[j][1]);
			}
		}
	}
}

void GA::printval() {
	//for (int i = 0;i < Vars;i++) {
	//	cout << population[N].gene[i] << " ";
	//}
	//cout << population[N].fitness << endl;
	avg[iii++] = population[N].fitness;
}
int main()
{
	srand((unsigned int)time(NULL));
	int ten = 10;
	while (ten--) {
		int T = MaxGens;
		GA ga;
		ga.initialize();
		while (T--) {
			ga.countfit();
			ga.countmaxfit();
			ga.RWS();
			ga.crossover();
			ga.variation();
		}
		ga.printval();
	}
	iii = 0;
	double aaa = 0.0;
	for (int i = 0;i < 10;i++) {
		aaa += avg[i];
	}
	aaa = aaa / 10.0;
	cout << "Pc= " << Pc << endl;
	cout << "Pm= " << Pm << endl;
	cout <<"average= "<< aaa << endl;

	
	for (double i = 0.1;i < 0.9;i = i + 0.1) {
		Pc = i;
		cout << "Pc= " << Pc << endl;
		for (double j = 0.01;j <= 0.09;j = j + 0.01) {
			Pm = j;
			int ten = 10;
			while (ten--) {
				int T = MaxGens;
				GA ga;
				ga.initialize();
				while (T--) {
					ga.countfit();
					ga.countmaxfit();
					ga.RWS();
					ga.crossover();
					ga.variation();
				}
				ga.printval();
			}
			iii = 0;
			double aaa = 0.0;
			for (int i = 0;i < 10;i++) {
				aaa += avg[i];
			}
			aaa = aaa / 10.0;
			cout << "Pm= " << Pm << ";  ";
			cout << "average= " << aaa << endl;
		}
	}
	return 0;
}

