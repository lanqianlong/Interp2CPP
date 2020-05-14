void locate(double xx[], int n, double x, int& j)
{
	int jl,ju,jm;
    jl = 0;
    ju = n + 1;
loop:  if (ju - jl > 1)
		{
        jm = (ju + jl) / 2;
        if ((xx[n] > xx[1]) && (x > xx[jm] ))
            jl = jm;
        else
		{
            ju = jm;
		}
        goto loop;
		}
    j = jl;
}

