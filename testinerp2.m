clc;
clear all;
%TestInterp2
[x,y]=meshgrid(0:2,0:2);
vl = [1:3;4:6;7:9];
[x2,y2]=meshgrid(1:0.1:2,0:2);
zz = interp2(x,y,vl,x2,y2,'spline');
% [x3,y3]=meshgrid(1:0.1:2,0:0.1:2);
% zz2 = interp2(x,y,vl,x3,y3,'bicubic');

%TestInterp2OnePoint
years = 1950:10:1990;
service = 10:10:30;
wage = [150.697 199.592 187.625
          179.323 195.072 250.287
          203.212 179.092 322.767
          226.505 153.706 426.730
          249.633 120.281 598.243];
w = interp2(service,years,wage,15,1975) 


%TestInterp2_extra
[X,Y] = meshgrid(-3:.25:3);
Z = peaks(X,Y);
[XI,YI] = meshgrid(-3:.125:3);
ZI = interp2(X,Y,Z,XI,YI);
mesh(X,Y,Z), hold, mesh(XI,YI,ZI+15)
hold off
axis([-3 3 -3 3 -5 20])