theta=linspace(0,2*pi,50);
theta = [theta, theta(1)];
x = 1/2*cos(theta);
y = 1/2*sin(theta)+1;

px = -2.5:0.1:2.5;
E = 1;
py = sqrt(E+px.^2);

hold off
plot(px,py,'r');
hold on
E = 1.5^2;
py = sqrt(E+px.^2);
plot(px,py,'r');
E = 0.5^2;
py = sqrt(E+px.^2);
plot(px,py,'r');


plot(x,y);
for t = 0.5:0.3:1.1
  qE = sqrt(x.^2+y.^2);
  C = asinh(x./qE);
  
  xx = qE.*sinh(t+C);
  yy = qE.*cosh(t+C);
  plot(xx,yy);
endfor

hold off
