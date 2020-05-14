void splin2(double x1a[], double x2a[], double ya[][11], double y2a[][11],
			int m, int n, double x1, double x2, double& y)
{ 
	double ytmp[100], y2tmp[100], yytmp[100];
	int j,k;
    for (j = 1; j<=m; j++)
	{
        for (k = 1; k<=n; k++)
		{
            ytmp[k] = ya[j][k];
            y2tmp[k] = y2a[j][k];
		}
        splint(x2a, ytmp, y2tmp, n, x2, yytmp[j]);
    }
    spline(x1a, yytmp, m, 1e+30, 1e+30, y2tmp);
    splint(x1a, yytmp, y2tmp, m, x1, y);
}
