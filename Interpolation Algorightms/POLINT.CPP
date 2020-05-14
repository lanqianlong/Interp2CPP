void polint(double xa[], double ya[], int n, double& x, double& y, double& dy)
{
    double c[11], d[11],dift,dif,ho,hp,w,den;
    int i,m,ns = 1;
    dif = fabs(x - xa[1]);
    for ( i = 1; i<=n; i++)
	{
        dift = fabs(x - xa[i]);
        if (dift < dif)
		{
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    y = ya[ns];
    ns = ns - 1;
    for(m = 1; m<=n-1; m++)
	{
        for (i = 1; i<=n-m; i++)
		{
            ho = xa[i] - x;
            hp = xa[i + m] - x;
            w = c[i + 1] - d[i];
            den = ho - hp;
            if (den == 0.0) 
			{
                cout<<"pause"<<endl;
                return;
			}
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
		}
		if (2 * ns < n - m )
            dy = c[ns + 1];
        else
		{
            dy = d[ns];
            ns = ns - 1;
		}
       
        y = y + dy;
	}
}
