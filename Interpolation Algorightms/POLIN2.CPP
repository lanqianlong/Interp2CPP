void polin2(double x1a[],double x2a[],double ya[][6],int m,int n,
			double x1, double x2,double&y, double& dy)
{
	double yntmp[20], ymtmp[20];
	int j,k;
    for(j = 1; j<=m; j++)
	{
        for(k = 1; k<=n; k++)
            yntmp[k] = ya[j][k];
        polint(x2a, yntmp, n, x2, ymtmp[j], dy);
    }
    polint(x1a, ymtmp, m, x1, y, dy);
}
