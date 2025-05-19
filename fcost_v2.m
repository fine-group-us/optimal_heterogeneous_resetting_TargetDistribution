function [F1,F2,G1,G2,cost]=fcost_v2(r,pT)
N=length(r)-1;
dx=2/N;
dx2=dx^2;
F1=zeros(1,N+1);
G1=F1;
F2=F1;
G2=F1;
G2(1)=1;
F1(2)=-0.5*dx^2;
F2(2)=dx;
for k=2:N
    F1(k+1)=-F1(k-1)+(2+dx2*r(k))*F1(k)-dx2;
    G1(k)=0.5*(F1(k+1)-F1(k-1))/dx;
    F2(k+1)=-F2(k-1)+(2+dx2*r(k))*F2(k);
    G2(k)=0.5*(F2(k+1)-F2(k-1))/dx;
end
F1r=-F1(N)+(2+dx2*r(N+1))*F1(N+1)-dx2;
G1(N+1)=0.5*(F1r-F1(N))/dx;
F2r=-F2(N)+(2+dx2*r(N+1))*F2(N+1);
G2(N+1)=0.5*(F2r-F2(N))/dx;

tau0=F1(1+N/2:N+1)-F1(1+N/2)/F2(1+N/2)*F2(1+N/2:N+1);
cost=-2*((tau0(1)*pT(1+N/2)+tau0(N/2+1)*pT(N+1))/2+tau0(2:N/2)*pT(2+N/2:N)')*dx;
end