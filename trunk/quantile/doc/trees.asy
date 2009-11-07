size(200,200);
real u=50;
real d=-50;
real l=100;

draw((0,0)--(l,u),Arrow);
draw((0,0)--(l,d),Arrow);
label("$x_0$",(0,0),W);
label("$x_0 + \mu$",(l,u),E);
label("$x_0 - \mu$",(l,d),E);
label("$p$",.5*(l,u),NW);
label("$1-p$", .5*(l,d),SW);
