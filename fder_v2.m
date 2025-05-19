function [mu]=fder_v2(r,pT,F1,F2,G1,G2)
N=length(r)-1;
dx=2/N;
dx2=dx^2;

F3=zeros(1,N+1);
G3=F3;
F4=F3;
G4=F3;
F3(2+N/2)=-0.5*dx2;
F3(N/2)=-0.5*dx2;
G4(1+N/2)=1;
F4(2+N/2)=dx;
F4(N/2)=-dx;

U12=zeros(N+1,N+1);
M=zeros(N+1,N+1);
for k=2+N/2:N
    %Computation of F4 from 0 to 1
    F4(k+1)=-F4(k-1)+(2+dx2*r(k))*F4(k);
    G4(k)=0.5*(F4(k+1)-F4(k-1))/dx;
    
    %Computation of F3,F4 from 0 to -1
    F3(1+N-k)=-F3(3+N-k)+(2+dx2*r(2+N-k))*F3(2+N-k)-dx2;
    G3(2+N-k)=0.5*(F3(3+N-k)-F3(1+N-k))/dx;
    F4(1+N-k)=-F4(3+N-k)+(2+dx2*r(2+N-k))*F4(2+N-k);
    G4(2+N-k)=0.5*(F4(3+N-k)-F4(1+N-k))/dx;
    %Computation of the propagators U12 and the functional derivative for x>0
    U12(k,k-1)=dx;
    M(k,k-1)=-U12(k,k-1)*(F1(k-1)-F1(1+N/2)*F2(k-1)/F2(1+N/2));
    for k1=k:N
        U12(k1+1,k-1)=-U12(k1-1,k-1)+(2+dx2*r(k1))*U12(k1,k-1);
        M(k1+1,k-1)=-U12(k1+1,k-1)*(F1(k-1)-F1(1+N/2)*F2(k-1)/F2(1+N/2));
    end
end
F4r=-F4(N)+(2+dx2*r(N+1))*F4(N+1);
G4(N+1)=0.5*(F4r-F4(N))/dx;
F3l=-F3(2)+(2+dx2*r(1))*F3(1)-dx2;
G3(1)=0.5*(F3(2)-F3l)/dx;
F4l=-F4(2)+(2+dx2*r(1))*F4(1);
G4(1)=0.5*(F4(2)-F4l)/dx;
U12(N+1,N)=dx;
M(N+1,N)=-U12(N+1,N)*(F1(N)-F1(1+N/2)*F2(N)/F2(1+N/2));

%Computation of the functional derivative for x<0
M(N/2+1:N+1,1:N/2+1)=-(F2(1:N/2+1).*F4(1:N/2+1).*(G1(1:N/2+1)-G3(1:N/2+1))+F2(1:N/2+1).*F3(1:N/2+1).*G4(1:N/2+1)-F1(1:N/2+1).*F4(1:N/2+1).*G2(1:N/2+1)).*F2(1:N/2+1)./(F4(1:N/2+1).*G2(1:N/2+1)-F2(1:N/2+1).*G4(1:N/2+1)).^2.*F4(N/2+1:N+1)';

%Symmetrization
M(1:N/2,1:N+1)=M(N+1:-1:N/2+2,N+1:-1:1);
mu=dx*((pT(1)*M(1,:)+pT(N+1)*M(N+1,:))/2+pT(2:N)*M(2:N,:));
end