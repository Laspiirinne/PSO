#ifndef _PSO_H_
#define _PSO_H_

#include<cmath>

/*各种适应度函数选择，要用哪个，就设置为1,但只能有一个为1*/
#define FUNC_CHOICE "Generalized_Penalized_2"                       //改！！！！！！！！！！
//二维函数
#define six_hump_camel_back_function 0
//单峰函数
#define Sphere 0							//range = 100.0,	Best_fitness = 0
#define Sphere_Shifted 0//range = 100.0,	Best_fitness = -450.0,	fbias = -450.0
#define Axis_parallel_hyper_ellipsoid 0
//多峰函数
#define Rastrigin	0				//range = 5.12,		Best_fitness = 0
#define Rastrigin_Shifted	0				//range = 5.0,		Best_fitness = -330.0,	fbias = -330.0 
#define Rastrigin_Shifted_Rotated	0		//range = 5.0,		Best_fitness = -330.0,	fbias = -330.0 
#define Ackley 0							//range = 32.0,		Best_fitness = 0
#define Ackley_Shifted_Rotated 0			//range = 32.0,		Best_fitness = -140.0,	fbias = -140.0
#define Griewank 0							//range = 600.0,	Best_fitness = 0
#define Schwefel	0  		//range = 500.0,	Best_fitness = 0
#define Generalized_Penalized_1 0			//range = 50.0,		Best_fitness = 0
#define Generalized_Penalized_2 1			//range = 50.0,		Best_fitness = 0
//非正式测试函数
#define Martin_and_Gaddy	0
#define Goldstein_and_Price 0
#define Demo 0
#define Rosenbrock 0						//range = 100.0,	Best_fitness = 0

#define PI acos(-1)
#define leftrange  -50.0                               //  改！！！！！！！
#define rightrange 50.0                            // 改！！！！！
#define Dim 30//粒子维度
#define PNum 40	//种群规模
#define ITE_N  10000//最大迭代次数
#define REALNUM 0.0    //函数的真实值                        改！！！！！！！
#define N 1//测试组数
#define L 50//每组测试数目
#define sudu 10.0    //最大速度                              改!!!!!!!!leftrange的1/5
extern int cur_n;			//当前迭代次数

							/*惯性权重函数*/
#define W_START 1.4
#define W_END	0.4
//#define W_FUN	(W_START-(W_START-W_END)*pow((double)cur_n/ITE_N,2))//0.9-cur_n/ITE_N*0.5
#define W_FUN 0.9-cur_n/ITE_N*0.5
/*个体和种群结构体*/
class PARTICLE
{
public:
	double X[Dim];
	double P[Dim];
	double V[Dim];
	double Fitness;
	friend class SWARM;
};
class SWARM
{
public:
	PARTICLE Particle[PNum];
	int GBestIndex;//全局最优粒子的编号
	double GBest[Dim];
	double Xup[Dim];
	double C1;
	double C2;
	double Xdown[Dim];
	double Vmax[Dim];
public:
	void RandInitofSwarm(void);
	void UpdateofVandX(void);
	void UpdatePandGbest(void);
};
double ComputAFitness(double X[]);


#endif


