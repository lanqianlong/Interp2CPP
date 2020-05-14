/*
Author:Qianlong Lan
*/
#pragma once
#ifndef INTERP2_H
#define INTERP2_H
#endif 
#include"math.h"
#include<iostream>
#include <string.h> 
#include<vector>
#ifndef SafeDeleteVec
#define SafeDeleteVec(X) { if((X)) delete (X); (X)=NULL;}
#endif
using namespace std;
class Interp
{
public:


	const char* strMethod1;
	const char* strMethod2;
	const char* strMethod3;


public:
	Interp();
	~Interp();
	/*****two interpolatino function**/
	void interp_onePoint(double* x, double* y, int n, double t, double* fval, const char* method);//return single point value
	void interp_multiPoint(double* x, double* y, int n, double* t, double* fval, int m, const char* method);//return multipoint value, point number is m
	double interp2_onePoint(double* x, double* y, double** z, int m, int n, double a, double b, const char* method);//matlab interp2 
	void interp2_multiPoint(double* x, double* y, double** z, int m, int n, double** a, double** b, int asize, int bsize, vector<vector<double>>&, const char* method);//matlab interp2 


private:
	//Lagrange linear interpolation
	double lgr(double* x, double* y, int n, double t);
	//Lagrange univariate three - point linear interpolation
	double lg3(double* x, double* y, int n, double t);
	//Smooth interpolation Cubic polynomial interpolation
	void spl(double* x, double* y, int n, int k, double t, double* s);

	int findIndex(double* x, double t, int n);
};




