void splie2(double x1a[], double x2a[], double ya[][11], int m, int n, double y2a[][11])
{
    double ytmp[101], y2tmp[101];
	int j,k;
    for (j = 1; j<=m; j++)
	{
        for (k = 1; k<=n; k++)
            ytmp[k] = ya[j][k];
        spline(x2a, ytmp, n, 1e+30, 1e+30, y2tmp);
        for (k = 1; k<=n; k++)
            y2a[j][k] = y2tmp[k];
    }
}
