#include"interp2.h"
#include<string.h>
Interp::Interp()
{
	strMethod1 = "lgr";
	strMethod2 = "lg3";
	strMethod3 = "spl";
}


Interp::~Interp()
{

}




double Interp::lgr(double* x, double* y, int n, double t)
{
	int i, j, k, m;
	double z, s;
	z = 0.0f;
	if (n < 1)
		return(z);
	if (n == 1)
	{
		z = y[0];
		return(z);
	}


	if (n == 2)
	{
		z = (y[0] * (t - x[1]) - y[1] * (t - x[0])) / (x[0] - x[1]);
		return(z);


	}
	if (n > 1 && t > x[n - 1]) {
		z = y[n - 1];
		return(z);
	}
	if (n > 1 && t < x[0]) {
		z = y[0];
		return(z);
	}




	i = 0;
	while ((x[i] < t) && (i < n)) {
		i = i + 1;
	}


	k = i - 4;
	if (k < 0)
		k = 0;
	m = i + 3;
	if (m > n - 1)
		m = n - 1;
	for (i = k; i <= m; i++)
	{
		s = 1.0;
		for (j = k; j <= m; j++) {
			if (j != i)
				s = s * (t - x[j]) / (x[i] - x[j]);
		}
		z = z + s * y[i];
	}
	return(z);


}



double Interp::lg3(double* x, double* y, int n, double t)
{
	int i, j, k, m;
	double z, s;
	z = 0.0;
	if (n < 1) return (z);
	if (n == 1) { z = y[0]; return(z); }
	if (n == 2)
	{
		z = (y[0] * (t - x[1]) - y[1] * (t - x[0])) / (x[0] - x[1]);
		return(z);
	}

	if (n > 0 && t < x[0]) {
		z = y[0];
		return(z);
	}

	if (n > 0 && t > x[n - 1]) {
		z = y[n - 1];
		return(z);
	}



	if (t <= x[1])
	{
		k = 0;
		m = 2;
	}
	else if (t >= x[n - 2])
	{
		k = n - 3;
		m = n - 1;
	}
	else
	{
		k = 1;
		m = n;
		while ((m - k) != 1)
		{
			i = (k + m) / 2;
			if (t < x[i - 1])
				m = i;
			else
				k = i;
		}
		k = k - 1;
		m = m - 1;
		if (fabs(t - x[k]) < fabs(t - x[m]))
			k = k - 1;
		else
			m = m + 1;
	}
	z = 0.0;
	for (i = k; i <= m; i++)
	{
		s = 1.0;
		for (j = k; j <= m; j++)
		{
			if (j != i)
				s = s * (t - x[j]) / (x[i] - x[j]);
		}
		z = z + s * y[i];




	}
	return (z);




}



void Interp::spl(double* x, double* y, int n, int k, double t, double* s)
{
	int kk, m, lc;


	double u[5], p, q;

	/*s[4] = 0.0f;
	s[0] = 0.0f;
	s[1] = 0.0f;
	s[2] = 0.0f;
	s[3] = 0.0f;*/

	s[4] = 0.0;
	s[0] = 0.0;
	s[1] = 0.0;
	s[2] = 0.0;
	s[3] = 0.0;




	if (n < 1)
		return;
	if (n == 1)
	{
		s[0] = y[0];
		s[4] = y[0];
		return;
	}
	if (n == 2)
	{
		s[0] = y[0];
		s[1] = (y[1] - y[0]) / (x[1] - x[0]);
		if (k < 0)
			s[4] = (y[0] * (t - x[1]) - y[1] * (t - x[0])) / (x[0] - x[1]);
		return;
	}


	if (k < 0 && n>0 && t < x[0])
	{


		s[4] = y[0];
		return;
	}


	if (k < 0 && n>0 && t > x[n - 1])
	{


		s[4] = y[n - 1];
		return;
	}






	if (k < 0)
	{


		if (t <= x[1])
			kk = 0;
		else if (t >= x[n - 1])
			kk = n - 2;
		else
		{


			kk = 1;
			m = n;
			while (((kk - m) != 1) && ((kk - m) != -1))
			{
				lc = (kk + m) / 2;
				if (t < x[lc - 1])
					m = lc;
				else
					kk = lc;


			}
			kk = kk - 1;
		}
	}
	else
		kk = k;


	if (kk > n - 1)
		kk = n - 2;


	u[2] = (y[kk + 1] - y[kk]) / (x[kk + 1] - x[kk]);
	if (n == 3)
	{
		if (kk == 0)
		{
			u[3] = (y[2] - y[1]) / (x[2] - x[1]);
			u[4] = 2.0 * u[3] - u[2];
			u[1] = 2.0 * u[2] - u[3];
			u[0] = 2.0 * u[1] - u[2];
		}
		else
		{
			u[1] = (y[1] - y[0]) / (x[1] - x[0]);
			u[0] = 2.0 * u[1] - u[2];
			u[3] = 2.0 * u[2] - u[1];
			u[4] = 2.0 * u[3] - u[2];
		}


	}
	else
	{
		if (kk <= 1)
		{
			u[3] = (y[kk + 2] - y[kk + 1]) / (x[kk + 2] - x[kk + 1]);
			if (kk == 1)
			{
				u[1] = (y[1] - y[0]) / (x[1] - x[0]);
				u[0] = 2.0 * u[1] - u[2];
				if (n == 4)
					u[4] = 2.0 * u[3] - u[2];
				else
					u[4] = (y[4] - y[3]) / (x[4] - x[3]);
			}
			else
			{
				u[1] = 2.0 * u[2] - u[3];
				u[0] = 2.0 * u[1] - u[2];
				u[4] = (y[3] - y[2]) / (x[3] - x[2]);
			}


		}
		else if (kk >= (n - 3))
		{
			u[1] = (y[kk] - y[kk - 1]) / (x[kk] - x[kk - 1]);
			if (kk == (n - 3))
			{
				u[3] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
				u[4] = 2.0 * u[3] - u[2];
				if (n == 4)
					u[0] = 2.0 * u[1] - u[2];
				else
					u[0] = (y[kk - 1] - y[kk - 2]) / (x[kk - 1] - x[kk - 2]);


			}
			else
			{
				u[3] = 2.0 * u[2] - u[1];
				u[4] = 2.0 * u[3] - u[2];
				u[0] = (y[kk - 1] - y[kk - 2]) / (x[kk - 1] - x[kk - 2]);
			}


		}
		else
		{
			u[1] = (y[kk] - y[kk - 1]) / (x[kk] - x[kk - 1]);
			u[0] = (y[kk - 1] - y[kk - 2]) / (x[kk - 1] - x[kk - 2]);
			u[3] = (y[kk + 2] - y[kk + 1]) / (x[kk + 2] - x[kk + 1]);
			u[4] = (y[kk + 3] - y[kk + 2]) / (x[kk + 3] - x[kk + 2]);
		}


	}


	s[0] = fabs(u[3] - u[2]);
	s[1] = fabs(u[0] - u[1]);
	if ((s[0] + 1.0 == 1.0) && (s[1] + 1.0 == 1.0))
		p = (u[1] + u[2]) / 2.0;
	else
		p = (s[0] * u[1] + s[1] * u[2]) / (s[0] + s[1]);
	s[0] = fabs(u[3] - u[4]);
	s[1] = fabs(u[2] - u[1]);
	if ((s[0] + 1.0 == 1.0) && (s[1] + 1.0 == 1.0))
		q = (u[2] + u[3]) / 2.0;
	else
		q = (s[0] * u[2] + s[1] * u[3]) / (s[0] + s[1]);


	s[0] = y[kk];
	s[1] = p;
	s[3] = x[kk + 1] - x[kk];
	s[2] = (3.0 * u[2] - 2.0 * p - q) / s[3];
	s[3] = (p + q - 2.0 * u[2]) / (s[3] * s[3]);
	if (k < 0)
	{
		p = t - x[kk];
		s[4] = s[0] + s[1] * p + s[2] * p * p + s[3] * p * p * p;

	}
	return;


}


double Interp::interp2_onePoint(double* x, double* y, double** z, int m, int n, double a, double b, const char* method)
{


	//find a，b
	int tempi, tempj;
	double w1, w2, w;
	//double* tempx(2), tempz(2);
	double tempx[2] = { 0,0 };
	double tempz[2] = { 0,0 };

	tempi = findIndex(x, a, m);
	tempj = findIndex(y, b, n);
	if (tempi < 0)
		tempi = 0;
	if (tempi > m - 2)
		tempi = m - 2;
	if (tempj < 0)
		tempj = 0;
	if (tempj > n - 2)
		tempj = n - 2;


	tempx[0] = x[tempi];
	tempx[1] = x[tempi + 1];
	/*tempz[0] = z[tempj * m + tempi];
	tempz[1] = z[tempj * m + tempi + 1];*/
	tempz[0] = z[tempj][tempi];
	tempz[1] = z[tempj][tempi + 1];

	//tempz[0] = z[tempi][tempj];
	//tempz[1] = z[tempi + 1][tempj];

	interp_onePoint(tempx, tempz, 2, a, &w1, method);
	//update tempz
	/*tempz[0] = z[(tempj + 1) * m + tempi];
	tempz[1] = z[(tempj + 1) * m + tempi + 1];*/
	tempz[0] = z[(tempj + 1)][tempi];
	tempz[1] = z[(tempj + 1)][tempi + 1];

	//tempz[0] = z[tempi][tempj+1];
	//tempz[1] = z[tempi+1][tempj + 1];
	interp_onePoint(tempx, tempz, 2, a, &w2, method);


	//update tempx and tempz
	tempx[0] = y[tempj];
	tempx[1] = y[tempj + 1];


	tempz[0] = w1;
	tempz[1] = w2;


	interp_onePoint(tempx, tempz, 2, b, &w, method);
	//}
	return(w);
}




void Interp::interp2_multiPoint(double* x, double* y, double** z, int m, int n, double** a, double** b, int asize, int bsize, vector<vector<double>>& fval, const char* method)
{
	/*
		Input x, y, out:z
		Input new a, b, out fval
	*/
	//grid nets
	for (int i = 0; i < asize; i++)
	{
		/*double test = interp2_onePoint(x, y, z, m, n, a[i], b[i], method);
		*(fval + i) = interp2_onePoint(x, y, z, m, n, a[i], b[i], method);*/
		for (int j = 0; j < bsize; j++) {
			double test = interp2_onePoint(x, y, z, m, n, a[i][j], b[i][j], method);

			fval[i][j] = interp2_onePoint(x, y, z, m, n, a[i][j], b[i][j], method);
			//fval.at(i).at(j) = interp2_onePoint(x, y, z, m, n, a[i], b[j], method);
		}

	}

}


/*
void Interp::interp2_multiPoint(double *x,double *y,double *z,int m,int n,double *gridx,double *gridy,int asize,int bsize,double *fval,char *method)
{
//grid nets
for(int i=0;i<asize;i++)
{
for(int j = 0;j<bsize;j++)
*(fval+j*asize+i) = interp2_onePoint(x,y,z,m,n,gridx[j*asize + i],gridx[j*asize + i],method);
}

}*/


int Interp::findIndex(double* x, double t, int n)
{
	int kk, m, lc;
	/*if(t<x[0])
	return( 0);
	if(t>=x[n-1])
	return(n-1);*/


	if (t <= x[1])
		kk = 0;
	else if (t >= x[n - 1])
		kk = n - 2;
	else
	{


		kk = 1;
		m = n;
		while (((kk - m) != 1) && ((kk - m) != -1))
		{
			lc = (kk + m) / 2;
			if (t < x[lc - 1])
				m = lc;
			else
				kk = lc;


		}
		kk = kk - 1;
	}


	return(kk);


}


void Interp::interp_onePoint(double* x, double* y, int n, double t, double* fval, const char* method)
{

	if (nullptr == method)
		*fval = lgr(x, y, n, t);
	else if (method == strMethod1)
		*fval = lgr(x, y, n, t);
	else if (method == strMethod2)
		*fval = lg3(x, y, n, t);
	else if (method == strMethod3)
	{
		static double s[5];
		spl(x, y, n, -1, t, s);
		*fval = s[4];
	}
	else
		*fval = lgr(x, y, n, t);
}


void Interp::interp_multiPoint(double* x, double* y, int n, double* t, double* fval, int m, const char* method)
{
	if (t == NULL)
		return;
	double tempVal = 0.0;;
	for (int k = 0; k < m; k++) {

		interp_onePoint(x, y, n, t[k], &tempVal, method);
		fval[k] = tempVal;
	}



}



