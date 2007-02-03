E = -1.2;
x= linspace(-pi*5,pi*7,200);
M = 0.15;

hold off
y = sqrt(E+cos(x)+M*x);
plot(x,y,'r',x,-y,'r');
hold on
for E = -1:0.4:1.6
  y = sqrt(2*(E+cos(x)+M*x));
  plot(x,y,'r',x,-y,'r');
endfor
hold off
