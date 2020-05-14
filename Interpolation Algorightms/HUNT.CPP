void hunt(double xx[], int n,double x, int& jlo)
{
	int lhi,inc,jm;
    int ascnd = xx[n] > xx[1];
    if( (jlo<= 0 )|| (jlo > n)) 
	{
        jlo = 0;
        lhi = n + 1;
        goto loop3;
	}
    inc= 1;
    if( (x >= xx[jlo] )&& ascnd )
	{
		loop1:   lhi = jlo + inc;
		if (lhi > n )
			lhi= n + 1;
		else 
			if ((x >= xx[lhi] )&& ascnd )
			{
				jlo= lhi;
				inc = inc + inc;
				goto loop1;
			}
	}
    else
	{
		lhi = jlo;
		loop2:   jlo = lhi - inc;
		if (jlo < 1 )
			jlo = 0;
		else 
		{
			if((x < xx[jlo] )&& ascnd )
			{
				lhi = jlo;
				inc = inc + inc;
				goto loop2;
			}
		}
     
    }
    loop3:  if(( lhi - jlo) == 1 )
			    return;
   	jm = (lhi + jlo)/ 2;
	if( (x > xx[jm] )&& ascnd)
		jlo= jm;
	else
		lhi = jm;
	goto loop3;
}

