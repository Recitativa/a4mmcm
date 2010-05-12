size(200,200);
real u=70;
real shr=0.5;
real w=50;

real l;
real h;
real s;
int i=0;
int ind;
  
int L=4;

pair[] pp;

pp.push((0,0));

for(ind=1; ind < 2^(L-1);++ind) {
  h = u*shr^floor(log(ind)/log(2));
  pp.push(pp[ind-1]+(w,h));
  pp.push(pp[ind-1]+(w,-h));
}

for(ind=1;ind < 2^L;++ind) {
  dot(pp[ind-1]);
  label(string(ind),pp[ind-1],N);
}

for(ind=2;ind < 2^L;++ind) {
  draw(pp[floor(ind/2)-1]--pp[ind-1]);
}

