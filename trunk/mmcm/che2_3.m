theta=linspace(0,2*pi,25);
theta = [theta, theta(1)];
x = 1/2*cos(theta);
y = 1/2*sin(theta)+1;

tn = 20;
eta = cos(linspace(0,pi,tn));

xa = 0; ya = 1.5;
E = ya^2/2-cos(xa);
bx = acos(-E);
tx = bx*eta;
px = [tx,tx(tn:-1:1)];
py = [sqrt(2*(E+cos(tx))),-sqrt(2*(E+cos(tx(tn:-1:1))))];
hold off
plot(px,py,'r');
hold on

xa = 0; ya = 0.5;
E = ya^2/2-cos(xa);
bx = acos(-E);
tx = bx*eta;
px = [tx,tx(tn:-1:1)];
py = [sqrt(2*(E+cos(tx))),-sqrt(2*(E+cos(tx(tn:-1:1))))];
plot(px,py,'r');


plot(x,y,";t=0;");

tt = 0;
for i = 1:6;
  t = 3;
  xx = [];
  yy = [];
  for i=1:length(x)
    [tx] = lsode("ff",[x(i),y(i)],[0,t],0);
    xx = [xx tx(2,1)];
    yy = [yy tx(2,2)];
  endfor
  tt += t;
  tstr = sprintf(";t=%d;",tt);
  plot(xx,yy,tstr);
  x = xx; y = yy;
endfor

hold off


function dx = ff(x,t)
  dx = [x(2), -sin(x(1))];
endfunction
