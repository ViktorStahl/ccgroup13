function U = BSeuCallspread_RBFFD(S,T,r,sig1,sig2,rho)

% The code expects a column vector, so if S is in a row, fix this
%
%---------   input  ----------------
% S - spot price
% T - time to maturity
% r - risk free interest rate
% sig1 - volatility on 1st asset
% sig2 - volatility on 2nd asset
% rho - correlation
%
%---------  output  ----------------
% U - option value

par={r,sig1,sig2,rho};

S2min=0; S2max=400;
S1min=0; S1max=4*S2max;

K=0;

%% Grid
Nx=300;

x=linspace(S1min,S1max,4*Nx)/S2max;
dx=x(2)-x(1);
y=linspace(S2min,S2max,Nx)/S2max;

xvec=[]; yvec=[]; indup=[]; inddiag=[];
for ii=1:numel(x)-2*Nx
    inddiag=[inddiag numel(xvec)+1];
    xl=x(ii)*ones(1,round((4*Nx-ii)/3+1));
    yl=linspace(0.25*x(ii),y(end),round((4*Nx-ii)/3+1));
    xvec=[xvec xl];
    yvec=[yvec yl];
    indup=[indup numel(xvec)];
end

indri=(numel(xvec)-round((4*Nx-ii)/3+1)+1):numel(xvec);

N=numel(xvec);

ind=1:N;
indcf=1:round((4*Nx-1)/3+1); 
% indcf=1;

indff=[inddiag(2:end-1),indri];
% indff=[indup(2:end),inddiag(2:end)];
% indff=[indup(1:end),inddiag(2:end)];

% indnd=[];
indnd=indup(2:end-1);

indin=ind; indin([indff,indcf,indnd])=[];

dt=0.001;
t=T:-dt:0;
M=length(t);

%% Initial condition
xvec=xvec'; yvec=yvec';
u0=max(xvec-yvec,zeros(length(xvec),1));
u=u0;

%% RBF
ep=0.00001;
n=7;
s=[xvec yvec];

%% Weights
W=sparse(N,N);

indc=findKNearestNeighbors(s,s,n);
% internal points {
for jj=1:numel(indin)
    ii=indin(jj);
    sc=[xvec(indc(ii,:)),yvec(indc(ii,:))];
    se=[xvec(ii),yvec(ii)];
    
    wc=rbfga_weights_BS2('spread',ep,sc,se,par);
    W(ii,indc(ii,:))=wc;
end
% } internal points

% boundary points {
for jj=1:numel(indnd)
    ii=indnd(jj);
    sc=[xvec(indc(ii,:)),yvec(indc(ii,:))];
    se=[xvec(ii),yvec(ii)];
    
    
    wc=rbfga_weights_BS2('yy0',ep,sc,se,par);
    W(ii,indc(ii,:))=wc;
end
% } boundary points

I=speye(size(W));

% save(['weightstri',num2str(N),'n',num2str(n),'.mat'])

%% Integration
% BDF-1
A=I-dt*W;
[L1 U1]=lu(A);

u1=u;
b=u1;
b(indcf)=0;
b(indff)=max(xvec(indff)-yvec(indff),zeros(numel(indff),1));
% b(indff)=max(xvec(indff)-yvec(indff)*exp(-r*dt),zeros(numel(indff),1));

u=L1\b;
u=U1\u;
u=max(u,zeros(size(u)));

% BDF-2
A=I-(2/3)*dt*W;
[L1 U1]=lu(A);

tri = delaunay(xvec,yvec);
for ii=3:M
    u2=u1;
    u1=u;
    b=(4/3)*u1-(1/3)*u2;
    b(indff)=max(xvec(indff)-yvec(indff),zeros(numel(indff),1));
    b(indcf)=0;
    
    u=L1\b;
    u=U1\u;
    
    u=max(u,zeros(size(u)));
end

u=u*S2max;
xvec=xvec*S2max;
yvec=yvec*S2max;

U=transpose(griddata(xvec,yvec,u,S(:,1),S(:,2),'cubic'));
end