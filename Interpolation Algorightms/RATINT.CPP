void ratint(double xa[],double  ya[], int n,double x, double& y, double& dy)
{
	int i,ns,m;
	double tiny,w,h,hh,dd,t,c[10], d[10];
    tiny = 1e-25;
    ns = 1;
    hh = fabs(x - xa[1]);
    for( i = 1; i<=n; i++)
	{
        h = fabs(x - xa[i]);
        if (h ==0.0 )
			{
				y = ya[i];
				dy = 0.0;
				return;
			}
        else if( h < hh )
			{
				ns = i;
				hh = h;
			}
        c[i] = ya[i];
        d[i] = ya[i] + tiny;
	}
    y = ya[ns];
    ns = ns - 1;
    for( m = 1; m<=n-1; m++)
	{
        for(i = 1; i<=n-m; i++)
		{
            w = c[i + 1] - d[i];
            h = xa[i + m] - x;
            t = (xa[i] - x) * d[i] / h;
            dd = t - c[i + 1];
            if (dd == 0.0)
				{
					 cout<<"pause"<<endl;
					 return;
				}
            dd = w / dd;
            d[i] = c[i + 1] * dd;
            c[i] = t * dd;
        }
        if (2 * double(ns) <double( n - m) )
		{
            dy = c[ns + 1];
		}
        else 
			
		{
            dy = d[ns];
            ns = ns - 1;
		}
        
        y = y + dy;
    }
}
