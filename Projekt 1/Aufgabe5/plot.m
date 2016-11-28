h=1/(d+1);
tx=ty=linspace(0,1,d+2);
B=shift(shift(resize(A,d+2),1),1,2);
surf(tx,ty,B);