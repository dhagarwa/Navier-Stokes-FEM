function [J,r] = jacres(p,e,t,u)
% triangle corner nodes
i=t(1,:); j=t(2,:); k=t(3,:);
% find triangle midpoints
x=(p(1,i)+p(1,j)+p(1,k))/3;
y=(p(2,i)+p(2,j)+p(2,k))/3;
% evaluate u, a, a’, and f
uu=(u(i)+u(j)+u(k))/3;
aa=a(uu);
da=(a(uu+1.e-8)-a(uu))/1.e-8;
ff=f(x,y);
% assemble jacobian and residual
[Aa ,unused,b]=assema(p,t,aa',0,ff);
[Ada,unused] =assema(p,t,da',0,0);
J=diag(Ada*u)+Aa;
r=Aa*u-b;
% enforce B.C.
for i=1:size(e,2)
n=e(1,i);
J(n,:)=0;
J(n,n)=1; r(n)=0;
end
