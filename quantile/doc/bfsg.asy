import graph;

size(300,300,IgnoreAspect);

file in=input("bfsg.dat").line();

real[][] a=in.dimension(0,0);
a = transpose(a);

real[] x = 1/a[0];
real[] y = a[1];

draw(graph(x,y),red, marker(scale(2)*unitcircle,red));

xaxis("$\Delta t = T/N$",BottomTop,LeftTicks);
yaxis("price",LeftRight,RightTicks);

     
     
