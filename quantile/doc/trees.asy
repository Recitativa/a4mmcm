size(200,200);
real u=50;
real d=30;
real l=100;

pair Nod1 = (l,d+u), Nod2 = (l,d-u);
draw((0,0)--Nod1,Arrow);
draw((0,0)--Nod2,Arrow);
label("$x_0$",(0,0),W);
label("$x_0 + d+u$",Nod1,E);
label("$x_0 + d-u$",Nod2,E);
label("$\frac{1}{2}$",.5*Nod1,NW);
label("$\frac{1}{2}$", .5*Nod2,SW);

