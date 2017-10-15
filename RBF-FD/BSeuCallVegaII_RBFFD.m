function [Vega] = BSeuCallVegaII_RBFFD(S,K,T,r,sig)
dsig = 0.9e-4;

sigu=sig+0.5*dsig;
sigd=sig-0.5*dsig;

%%
N=12000;
N1=round(N/4);
N2=N-N1;

lim=4;
cen=0.97;

chebr=flip(chebx([1,0],N1)+2);
chebl=chebx([0,1],N2+1);

x1=chebr*cen;
x2=chebl*(lim-cen)+cen;

x=[x1;x2(2:end)];

dx=x(2)-x(1);
n=5; %stencil size
m=round((n-1)/2);

Ks=1;

M=2000;
dt=T/(M-1);
t=T:-dt:0;

% Initial condition
u=max(x-Ks,zeros(N,1)); 

% RBF
ep=0.001;

% Weights
W=sparse(N,N);

for ii=2:m
    indc=1:n;
    xc=x(ii);
    xi=x(indc);

    wc1=rbfga_weights_BS('xx',ep,xi,xc,r,sigu);
    wc1=wc1*0.5*sigu^2*x(ii)^2;
    
    wc2=rbfga_weights_BS('x',ep,xi,xc,r,sigu);
    wc2=wc2*r*x(ii);
    
    wc3=rbfga_weights_BS('0',ep,xi,xc,r,sigu);
    wc3=-wc3*r;
    
    wc=wc1+wc2+wc3;

    W(ii,indc)=wc;
end

for ii=m+1:N-m
    indc=ii-m:ii+m;
    xc=x(ii);
    xi=x(indc);

    wc1=rbfga_weights_BS('xx',ep,xi,xc,r,sigu);
    wc1=wc1*0.5*sigu^2*x(ii)^2;
    
    wc2=rbfga_weights_BS('x',ep,xi,xc,r,sigu);
    wc2=wc2*r*x(ii);
    
    wc3=rbfga_weights_BS('0',ep,xi,xc,r,sigu);
    wc3=-wc3*r;
    
    wc=wc1+wc2+wc3;

    W(ii,indc)=wc;
end

for ii=N-m+1:N-1
    indc=N-n+1:N;
    xc=x(ii);
    xi=x(indc);

    wc1=rbfga_weights_BS('xx',ep,xi,xc,r,sigu);
    wc1=wc1*0.5*sigu^2*x(ii)^2;
    
    wc2=rbfga_weights_BS('x',ep,xi,xc,r,sigu);
    wc2=wc2*r*x(ii);
    
    wc3=rbfga_weights_BS('0',ep,xi,xc,r,sigu);
    wc3=-wc3*r;
    
    wc=wc1+wc2+wc3;

    W(ii,indc)=wc;
end

% Integration
I=speye(N);

%BDF-1
A=I-W*dt;
[L,U]=lu(A);

P=2;

for ii=2:P    
    u1=u;
    
    b=u1;
    b(end)=x(end)-Ks*exp(-r*(ii-1)*dt);
    
    u=U\(L\b);
    u=max(u,0);
end

%BDF-2
A=I-(2/3)*dt*W;
[L,U]=lu(A);

for ii=P+1:M
    u2=u1;
    u1=u;
    
    b=((4/3)*u1-(1/3)*u2);
    b(end)=x(end)-Ks*exp(-r*(ii-1)*dt);
    
    u=U\(L\b);
    u=max(u,0);
end

x=K*x;
u=K*u;

U=interp1(x,u,S,'spline');

Uu=U;

%%
% N=8000;
N1=round(N/4);
N2=N-N1;

lim=4;
cen=0.97;

chebr=flip(chebx([1,0],N1)+2);
chebl=chebx([0,1],N2+1);

x1=chebr*cen;
x2=chebl*(lim-cen)+cen;

x=[x1;x2(2:end)];

dx=x(2)-x(1);
n=5; %stencil size
m=round((n-1)/2);

Ks=1;

% M=1500;
dt=T/(M-1);
t=T:-dt:0;

% Initial condition
u=max(x-Ks,zeros(N,1)); 

% RBF
ep=0.001;

% Weights
W=sparse(N,N);

for ii=2:m
    indc=1:n;
    xc=x(ii);
    xi=x(indc);

    wc1=rbfga_weights_BS('xx',ep,xi,xc,r,sigd);
    wc1=wc1*0.5*sigd^2*x(ii)^2;
    
    wc2=rbfga_weights_BS('x',ep,xi,xc,r,sigd);
    wc2=wc2*r*x(ii);
    
    wc3=rbfga_weights_BS('0',ep,xi,xc,r,sigd);
    wc3=-wc3*r;
    
    wc=wc1+wc2+wc3;

    W(ii,indc)=wc;
end

for ii=m+1:N-m
    indc=ii-m:ii+m;
    xc=x(ii);
    xi=x(indc);

    wc1=rbfga_weights_BS('xx',ep,xi,xc,r,sigd);
    wc1=wc1*0.5*sigd^2*x(ii)^2;
    
    wc2=rbfga_weights_BS('x',ep,xi,xc,r,sigd);
    wc2=wc2*r*x(ii);
    
    wc3=rbfga_weights_BS('0',ep,xi,xc,r,sigd);
    wc3=-wc3*r;
    
    wc=wc1+wc2+wc3;

    W(ii,indc)=wc;
end

for ii=N-m+1:N-1
    indc=N-n+1:N;
    xc=x(ii);
    xi=x(indc);

    wc1=rbfga_weights_BS('xx',ep,xi,xc,r,sigd);
    wc1=wc1*0.5*sigd^2*x(ii)^2;
    
    wc2=rbfga_weights_BS('x',ep,xi,xc,r,sigd);
    wc2=wc2*r*x(ii);
    
    wc3=rbfga_weights_BS('0',ep,xi,xc,r,sigd);
    wc3=-wc3*r;
    
    wc=wc1+wc2+wc3;

    W(ii,indc)=wc;
end

% Integration
I=speye(N);

%BDF-1
A=I-W*dt;
[L,U]=lu(A);

P=2;

for ii=2:P    
    u1=u;
    
    b=u1;
    b(end)=x(end)-Ks*exp(-r*(ii-1)*dt);
    
    u=U\(L\b);
    u=max(u,0);
end

%BDF-2
A=I-(2/3)*dt*W;
[L,U]=lu(A);

for ii=P+1:M
    u2=u1;
    u1=u;
    
    b=((4/3)*u1-(1/3)*u2);
    b(end)=x(end)-Ks*exp(-r*(ii-1)*dt);
    
    u=U\(L\b);
    u=max(u,0);
end

x=K*x;
u=K*u;

U=interp1(x,u,S,'spline');

Ul=U;

%%
Vega=(Uu-Ul)/dsig;
end