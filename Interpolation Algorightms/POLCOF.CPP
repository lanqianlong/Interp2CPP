void polcof(double xa[],double ya[],int n, double cof[])
{
	double x[15], y[15],xmin,dy,kk;
	int j,i,k;
	kk=0.0;
    for (j = 1; j<=n; j++)
	{
        x[j] = xa[j];
        y[j] = ya[j];
    }
    for (j = 1; j<=n; j++)
	{
        polint(x, y, n + 1 - j, kk, cof[j], dy);
        xmin = 1e+38;
        k = 0;
        for( i = 1; i<=n + 1 - j; i++)
		{
			
            if (fabs(x[i]) < xmin )
			{
                xmin = fabs(x[i]);
                k = i;
            }
            if (x[i] !=0.0)
				y[i] = (y[i] - cof[j]) / x[i];
        }
        if (k < n + 1 - j)
		{
            for (i = k + 1; i<=n + 1 - j; i++)
			{
                y[i - 1] = y[i];
                x[i - 1] = x[i];
            }
        }
    }
}
   
