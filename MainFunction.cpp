#pragma warning(disable:4996)
#include "PSO.h"
#include"stdlib.h"
#include"time.h"
#include "stdio.h"
#include"math.h"

int cur_n;
//struct PARTICLE particle;

int main(int argc, const char *argv[])
{
	double Maxfitness, Minfitness, Sumfitness, Avefitness, Maxcur_n, Mincur_n, Sumcur_n, Avecur_n, Sumtime;//基于整个文件
	FILE *fp1, *fp2;
	struct tm *newtime;
	char tmpbuffile1[128], tmpbuffile2[128], tmpbuf2file1[128], tmpbuf2file2[128];
	time_t lt1;
	lt1 = time(NULL);
	newtime = localtime(&lt1);                //得到现在时间
	int k, i; int endcur_n, Endcur_n;//终结时的迭代次数  
	int flag;               // 使阈值处只进入一次   标志是否到达阈值
	clock_t start, end, start3, end3;                          //start3和end3用来检验是否丢了一次运行的时间
	double sumfitness1, sumfitness2, sumtime, sumcur_n;      //sumtime是每组的运行时间
	double maxfitness1, maxfitness2, minfitness1, minfitness2;                 //最优适应度，最差适应度
	double maxcur_n, mincur_n;                     //最优迭代次数，最差迭代次数
	double avefitness1, avefitness2, avecur_n;					//平均适应度，平均迭代次数
	int bianhao, cishu;                           //为了满足阈值时继续循环        记录有多少次达到阀值
	double shiyingdu, wucha;
	srand((unsigned)time(NULL));
	sprintf(tmpbuffile1, "psoDIM_%dSNUM_%dITE_%d_%s%d%02d%02d_%d%02d%02d_file1.txt",
		Dim, PNum, ITE_N, FUNC_CHOICE,
		1900 + newtime->tm_year, 1 + newtime->tm_mon, newtime->tm_mday, newtime->tm_hour, newtime->tm_min, newtime->tm_sec);//输入文件的文件名
	sprintf(tmpbuffile2, "psoDIM_%dSNUM_%dITE_%d_%s%d%02d%02d_%d%02d%02d_file2.txt",
		Dim, PNum, ITE_N, FUNC_CHOICE,
		1900 + newtime->tm_year, 1 + newtime->tm_mon, newtime->tm_mday, newtime->tm_hour, newtime->tm_min, newtime->tm_sec);//输入文件的文件名
	sprintf(tmpbuf2file1, "psoDIM_%dSNUM_%dITE_%d_%s%d%02d%02d_%d%02d%02d_file1",
		Dim, PNum, ITE_N, FUNC_CHOICE,
		1900 + newtime->tm_year, 1 + newtime->tm_mon, newtime->tm_mday, newtime->tm_hour, newtime->tm_min, newtime->tm_sec);
	sprintf(tmpbuf2file2, "psoDIM_%dSNUM_%dITE_%d_%s%d%02d%02d_%d%02d%02d_file2",
		Dim, PNum, ITE_N, FUNC_CHOICE,
		1900 + newtime->tm_year, 1 + newtime->tm_mon, newtime->tm_mday, newtime->tm_hour, newtime->tm_min, newtime->tm_sec);
	fp1 = fopen(tmpbuffile1, "wt");
	fp2 = fopen(tmpbuffile2, "wt");
	fprintf(fp1, "%s\n\n\n", tmpbuf2file1);
	fprintf(fp2, "%s\n\n\n", tmpbuf2file2);
	fprintf(fp1, "粒子维度：%d\n种群规模:%d\n迭代次数:%d\n坐标范围左:%lf\n坐标范围右:%lf\n每组测试数目:%d\n测试组数:%d\n测试函数名称:%s\n\n\n", Dim, PNum, ITE_N, leftrange, rightrange, L, N, FUNC_CHOICE);   //在文件头部输入基本信息
	fprintf(fp2, "粒子维度：%d\n种群规模:%d\n迭代次数:%d\n坐标范围左:%lf\n坐标范围右:%lf\n每组测试数目:%d\n测试组数:%d\n测试函数名称:%s\n\n\n", Dim, PNum, ITE_N, leftrange, rightrange, L, N, FUNC_CHOICE);   //在文件头部输入基本信息
	Sumfitness = 0; Sumcur_n = 0; Sumtime = 0, Maxfitness = -9999999; Minfitness = 9999999; Maxcur_n = -9999999; Mincur_n = 9999999;
	for (k = 0; k < N; k++)
	{
		fprintf(fp1, "======================================================================================================================\n");
		fprintf(fp2, "======================================================================================================================\n");
		fprintf(fp1, "第%d组:\n\n", k + 1);
		fprintf(fp2, "第%d组:\n\n", k + 1);
		sumtime = 0; sumfitness1 = 0; sumfitness2 = 0; sumcur_n = 0; cishu = 0;
		maxfitness1 = -9999999; maxfitness2 = -9999999; minfitness1 = 9999999; minfitness2 = 9999999; maxcur_n = -9999999, mincur_n = 9999999;// start3 = clock();
		fprintf(fp2, "测试次数\t达到误差允许阀值时的代数\t粒子编号\t函数值\t函数误差\t迭代终止时的迭代次数\t粒子编号\t函数值\t函数误差\n");
		fprintf(fp2, "注：达到误差允许阀值时的代数 列中数据为 2000时为未在最大迭代次数限制内找到误差允许的值。\n");
		start = clock();
		for (i = 0; i < L; i++)
		{
			SWARM swarm;
			fprintf(fp1, "代数\t编号\t函数值\t函数误差\t\t坐标\n");
			cur_n = 0;
			flag = 1;

			swarm.RandInitofSwarm();             //初始化粒子群
			while (++cur_n <= ITE_N)
			{
				endcur_n = cur_n;
				swarm.UpdateofVandX(); //速度和位置更新，即飞翔
				swarm.UpdatePandGbest();    //更新个体历史最优解P和全局最优解GBest
				fprintf(fp1, "%d\t%d\t%E\t%E\t", cur_n, swarm.GBestIndex, ComputAFitness(swarm.GBest), ComputAFitness(swarm.GBest) - REALNUM);
				if (fabs(ComputAFitness(swarm.GBest) - REALNUM) < 1E-6&&flag == 1)
				{
					flag = 0;
					fprintf(fp1, "到达阀值");
					Endcur_n = endcur_n;
					bianhao = swarm.GBestIndex;                                    //带1的是真正意义上的最优解，带2的是程序运行完时的最优解
					shiyingdu = ComputAFitness(swarm.GBest);
					wucha = shiyingdu - REALNUM;
					cishu++;
				}
				/*for (int i = 0; i < Dim; i++)
				fprintf(fp1, "\t%E", swarm.GBest[i]);*/
				fprintf(fp1, "\n");

			}
			if (flag == 1)            //若未到达阀值，则进入
			{
				Endcur_n = endcur_n;
				bianhao = swarm.GBestIndex;
				shiyingdu = ComputAFitness(swarm.GBest);
				wucha = shiyingdu - REALNUM;
			}
			maxfitness1 = maxfitness1 > shiyingdu ? maxfitness1 : shiyingdu;
			minfitness1 = minfitness1 < shiyingdu ? minfitness1 : shiyingdu;
			maxfitness2 = maxfitness2 >  ComputAFitness(swarm.GBest) ? maxfitness2 : ComputAFitness(swarm.GBest);        //ComputAFitness(swarm.GBest)和swarm.particle[GBestIndex].fitness是不一样的
			minfitness2 = minfitness2 <  ComputAFitness(swarm.GBest) ? minfitness2 : ComputAFitness(swarm.GBest);
			maxcur_n = maxcur_n > Endcur_n ? maxcur_n : Endcur_n;
			mincur_n = mincur_n < Endcur_n ? mincur_n : Endcur_n;
			sumfitness1 += shiyingdu;
			sumfitness2 += ComputAFitness(swarm.GBest);
			sumcur_n += Endcur_n;
			if (Endcur_n < ITE_N)
			{
				fprintf(fp1, "第%-2d次测试:\n在\t%d\t代函数值达到误差允许：函数值:\t%E\t函数误差：\t%E\n", i + 1, Endcur_n, shiyingdu, wucha);
			}
			else
			{
				fprintf(fp1, "第%-2d次测试:\n在\t%d\t代停止未能在最大迭代次数内找到误差允许的函数值\t函数值:\t%E\t函数误差：\t%E\n", i + 1, Endcur_n, shiyingdu, wucha);
			}
			fprintf(fp1, "在\t%d\t代达到终止条件：领头羊的编号\t%d\t函数值：\t%E\t函数误差:\t%E\n\n", ITE_N, swarm.GBestIndex, ComputAFitness(swarm.GBest), ComputAFitness(swarm.GBest) - REALNUM);
			fprintf(fp2, "%-2d\t%-5d\t%-2d\t%E\t%E\t%-5d\t%d\t%E\t%E\n", i + 1, Endcur_n, bianhao, shiyingdu, wucha, ITE_N, swarm.GBestIndex, ComputAFitness(swarm.GBest), ComputAFitness(swarm.GBest) - REALNUM);
			//fprintf(fp, "第%-2d次测试:\n达到误差允许时的迭代次数\t%-5d\t粒子编号:\t%-2d\t函数值:\t%E\t函数误差:\t%E\n达到终止条件时的迭代次数:\t%-5d\t粒子编号:\t%d\t函数值：\t%E\t函数误差：\t%E\n\n", i + 1, Endcur_n, bianhao, shiyingdu,wucha,ITE_N, swarm.GBestIndex, ComputAFitness(swarm.GBest), ComputAFitness(swarm.GBest)-REALNUM);
			printf("组数：%d/%d  次数：%d/%dis OK!\n", k + 1, N, i + 1, L);
		}//end3 = clock();
		end = clock();
		avefitness1 = sumfitness1 / L;
		avefitness2 = sumfitness2 / L;
		avecur_n = sumcur_n / L;
		sumtime = (end - start) / CLOCKS_PER_SEC;
		if (cishu == 0)
		{
			fprintf(fp1, "未能在有限迭代次数内找到误差允许范围内的值。\n\n");
			fprintf(fp1, "每次测试都迭代到最大迭代次数的统计数据如下：\n");
			fprintf(fp1, "最优函数值，函数误差:\t%E\t%E\n", minfitness2, minfitness2 - REALNUM);
			fprintf(fp1, "最差函数值，函数误差:\t%E\t%E\n", maxfitness2, maxfitness2 - REALNUM);
			fprintf(fp1, "平均函数值，函数误差:\t%E\t%E\n", avefitness2, avefitness2 - REALNUM);
			fprintf(fp1, "运行时间:\t%.3f 秒\n\n", sumtime); //fprintf(fp, "%.3f", (double)(end3 - start3) / CLOCKS_PER_SEC);
			fprintf(fp1, "=======================================================================================================\n\n\n\n");
			fprintf(fp2, "\n未能在有限迭代次数内找到误差允许范围内的值。\n\n");
			fprintf(fp2, "每次测试都迭代到最大迭代次数的统计数据如下：\n");
			fprintf(fp2, "最优函数值，函数误差:\t%E\t%E\n", minfitness2, minfitness2 - REALNUM);
			fprintf(fp2, "最差函数值，函数误差:\t%E\t%E\n", maxfitness2, maxfitness2 - REALNUM);
			fprintf(fp2, "平均函数值，函数误差:\t%E\t%E\n", avefitness2, avefitness2 - REALNUM);
			fprintf(fp2, "运行时间:\t%.3f 秒\n\n", sumtime); //fprintf(fp, "%.3f", (double)(end3 - start3) / CLOCKS_PER_SEC);
			fprintf(fp2, "=======================================================================================================\n\n\n\n");
			printf("组数：%d/%dis OK\n", k + 1, N);
		}
		else
		{
			fprintf(fp1, "在\t%d\t次测试中，有\t%d\t次测试找到了误差允许的值，本组测试数据综合统计如下：\n", L, cishu);
			fprintf(fp1, "最优函数值，函数误差:\t%E\t%E\n", minfitness1, minfitness1 - REALNUM);
			fprintf(fp1, "最差函数值，函数误差:\t%E\t%E\n", maxfitness1, maxfitness1 - REALNUM);
			fprintf(fp1, "平均函数值，函数误差:\t%E\t%E\n", avefitness1, avefitness1 - REALNUM);
			fprintf(fp1, "最小迭代次数\t%.1f\n", mincur_n);
			fprintf(fp1, "最大迭代次数\t%.1f\n", maxcur_n);
			fprintf(fp1, "平均迭代次数\t%.1f\n", avecur_n);
			fprintf(fp1, "每次测试都迭代到最大迭代次数的统计数据如下：\n");
			fprintf(fp1, "最优函数值，函数误差:\t%E\t%E\n", minfitness2, minfitness2 - REALNUM);
			fprintf(fp1, "最差函数值，函数误差:\t%E\t%E\n", maxfitness2, maxfitness2 - REALNUM);
			fprintf(fp1, "平均函数值，函数误差:\t%E\t%E\n", avefitness2, avefitness2 - REALNUM);
			fprintf(fp1, "运行时间:\t%.3f 秒\n\n", sumtime); //fprintf(fp, "%.3f", (double)(end3 - start3) / CLOCKS_PER_SEC);
			fprintf(fp1, "=======================================================================================================\n\n\n\n");
			fprintf(fp2, "\n在\t%d\t次测试中，有\t%d\t次测试找到了误差允许的值，本组测试数据综合统计如下：\n", L, cishu);
			fprintf(fp2, "最优函数值，函数误差:\t%E\t%E\n", minfitness1, minfitness1 - REALNUM);
			fprintf(fp2, "最差函数值，函数误差:\t%E\t%E\n", maxfitness1, maxfitness1 - REALNUM);
			fprintf(fp2, "平均函数值，函数误差:\t%E\t%E\n", avefitness1, avefitness1 - REALNUM);
			fprintf(fp2, "最小迭代次数\t%.1f\n", mincur_n);
			fprintf(fp2, "最大迭代次数\t%.1f\n", maxcur_n);
			fprintf(fp2, "平均迭代次数\t%.1f\n", avecur_n);
			fprintf(fp2, "每次测试都迭代到最大迭代次数的统计数据如下：\n");
			fprintf(fp2, "最优函数值，函数误差:\t%E\t%E\n", minfitness2, minfitness2 - REALNUM);
			fprintf(fp2, "最差函数值，函数误差:\t%E\t%E\n", maxfitness2, maxfitness2 - REALNUM);
			fprintf(fp2, "平均函数值，函数误差:\t%E\t%E\n", avefitness2, avefitness2 - REALNUM);
			fprintf(fp2, "运行时间:\t%.3f 秒\n\n", sumtime); //fprintf(fp, "%.3f", (double)(end3 - start3) / CLOCKS_PER_SEC);
			fprintf(fp2, "=======================================================================================================\n\n\n\n");
			printf("组数：%d/%dis OK\n", k + 1, N);
		}
		Maxfitness = Maxfitness > maxfitness1 ? Maxfitness : maxfitness1;
		Minfitness = Minfitness < minfitness1 ? Minfitness : minfitness1;
		Maxcur_n = Maxcur_n > maxcur_n ? Maxcur_n : maxcur_n;
		Mincur_n = Mincur_n < mincur_n ? Mincur_n : mincur_n;
		Sumfitness += avefitness1;
		Sumcur_n += avecur_n;
		Sumtime += sumtime;
	}
	/*Avefitness = Sumfitness / N;
	Avecur_n = Sumcur_n / N;
	fprintf(fp, "\n\n\n\n\n全体最优适应度: %E\n", Minfitness);
	fprintf(fp, "全体最差适应度: %E\n", Maxfitness);
	fprintf(fp, "全体平均适应度: %E\n", Avefitness);
	fprintf(fp, "全体最优误差: %E\n", Minfitness-REALNUM);
	fprintf(fp, "全体最差误差: %E\n", Maxfitness-REALNUM);
	fprintf(fp, "全体平均误差: %E\n", Avefitness-REALNUM);
	fprintf(fp, "全体最优迭代次数：%.1f\n", Mincur_n);
	fprintf(fp, "全体最差迭代次数：%.1f\n", Maxcur_n);
	fprintf(fp, "全体平均迭代次数：%.1f\n", Avecur_n);
	fprintf(fp, "全体运行时间: %.3f 秒\n\n\n\n", Sumtime);
	fclose(fp);*/
	return 0;
}