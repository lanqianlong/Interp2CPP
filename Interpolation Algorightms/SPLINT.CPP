void splint(double xa[], double ya[], double y2a[], int n, double& x, double& y)
{
	int klo,khi,k;
	double h,a,b,aaa,bbb;
    klo = 1;
    khi = n;
loop:   if (khi - klo > 1 )
	{
		k = (khi + klo) / 2;
		if (xa[k] > x)
			khi = k;
		else
		{
			klo = k;
		}
           goto loop;
	}
    h = xa[khi] - xa[klo];
    if (h == 0 )
	{
		cout<<"  pause  'bad  xa  input'"<<endl;
		return;
	}
	a = (xa[khi] - x) / h;
	b = (x - xa[klo]) / h;
	aaa = a * ya[klo] + b * ya[khi];
	bbb = (a*a*a - a) * y2a[klo] + (b*b*b - b) * y2a[khi];
	y = aaa + bbb * h*h /6.0;
}
