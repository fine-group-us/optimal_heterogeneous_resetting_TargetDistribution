tic

%% Parameters and initialisation
npasos=1e6;
freq=100;
landa0=1e4;
N=2000;
dx=2/N;
x=-1:dx:1;

% Distribution that is used 
b=-0.75;
pT = (3-2*b)/6 + b*x.^2;

costv=zeros(npasos/freq+1,1);
muv=zeros(npasos/freq+1,N+1);
rv=zeros(npasos/freq+1,N+1);
landav=zeros(npasos/freq+1,1);

landa=0;
[F1,F2,G1,G2,cost]=fcost_v2(r,pT); % Obtain the average MFPT
[mu]=fder_v2(r,pT,F1,F2,G1,G2); % Obtain the functional derivative    
for cont=1:npasos+1
    
    % After freq steps, we save the current state 
    if mod(cont-1,freq)==0
    muv(1+(cont-1)/freq,:)=mu;    
    rv(1+(cont-1)/freq,:)=r;
    landav(1+(cont-1)/freq,:)=landa;
    costv(1+(cont-1)/freq)=cost;
    disp((cont-1)/freq)
    % figure(1)
    %plot(x(1+N/2:N+1),r(1+N/2:N+1))
    % figure(2)
    %plot(landav(landav>0),'o')
    % figure(3)
    %plot(costv(costv>0),'o')
    % figure(4)
    %plot(x(1+N/2:N+1),mu(1+N/2:N+1))
    end
    
    landa=landa0;
    costnew=cost+1;
    % If the new averaged MFPT is greater than the previous one, 
    % we reduce the step 
    while costnew>cost
        rnew=max(r-landa*mu,0); 
        [F1,F2,G1,G2,costnew]=fcost_v2(rnew,pT);
        landa=landa/2;
    end

    landa=landa*2;
    r=rnew;
    cost=costnew;   
    % Compute the new functional derivative (perturbation of r)
    [mu]=fder_v2(r,pT,F1,F2,G1,G2); 

end
toc