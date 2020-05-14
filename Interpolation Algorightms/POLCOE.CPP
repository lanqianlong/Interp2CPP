void polcoe(double x[], double y[], int n, double cof[])
{
    double phi,ff,b,s[15];
	int i,j,k;
    for(i = 1; i<=n; i++)
	{
		s[i] = 0.0;
		cof[i] = 0.0;
	}
    s[n] = -x[1];
    for (i = 2; i<=n; i++)
	{
		for (j = n + 1 - i;j<=n-1;j++)
		{
			s[j] = s[j] - x[i] * s[j + 1];
		}
		s[n] = s[n] - x[i];
	}
    for (j = 1; j<=n; j++)
	{
		phi = n;
		for (k = n - 1; k>=1; k--)
		{
			phi = k * s[k + 1] + x[j] * phi;
		}
		ff = y[j] / phi;
		b = 1.0;
		for (k = n; k>=1; k--)
		{
			cof[k] = cof[k] + b * ff;
			b = s[k] + x[j] * b;
		}
	}
}
