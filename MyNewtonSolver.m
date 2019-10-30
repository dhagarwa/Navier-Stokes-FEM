function MyNewtonSolver(geom)
[p,e,t]=initmesh(geom,'hmax',0.1);
u=zeros(size(p,2),1);
for k=1:5
[J,r]=jacres(p,e,t,u);
d=J\r;
u=u+d;
sprintf('|d|=%f, |r|=%f', norm(d), norm(r))
end
pdesurf(p,t,u)