function sigma = locvolII(ir,s0,ss,tt)
% This function evaluates a local volitility surface sigma(s,t) 
%  
% ir is the risk free interest rate (or could be carry cost)
%  
% s0 is a scalar value and should be the stock price that was used for
%    calculating the local volatility surface as well as the initial
%    stock value for which we are evaluating the price.  
%  
% ss is the stock price(s) for which we want to evaluate sigma.
%    ss can be a scalar, a vector, or a matrix.
%  
% tt is the time for which we want to evaluate sigma. If ss is a scalar,
%    tt can be a scalar, vector or matrix. If ss is a vector, tt can be a
%    scalar or a vector of the same size. If ss is a matrix, tt can be a
%    scalar or a matrix of the same size. sigma(s,t) will be computed for
%    the coordinate pairs indicated by the sizes of s and t.
%  
% NOTE: The function can not return a result for s=0, but gives
% reasonable results for rather small values of s and for times <=0.5.
  %
  % The global variance surface is given in an SVI parametrization
  % form. We have implemented one particular example here.
  % We start by generating functions for parameters and subexpressions 
  %
  % Coefficient functions and their derivatives
  %
  a = @(T) 0.01+0.03*sqrt(T+0.04);
  da = @(T) 0.015./sqrt(T+0.04);

  r = @(T) 0.06*(1-0.87*sqrt(T));
  dr = @(T) -0.06*0.87/2./sqrt(T);

  l = @(T) 0.31*(1-0.7*sqrt(T));
  dl = @(T) -0.31*0.7/2./sqrt(T);

  m = @(T) 0.03+0.01*T;
  dm = 0.01;

  q = @(T) 0.15*(0.4+0.6*sqrt(T+0.04));
  dq = @(T) 0.15*0.6/2./sqrt(T+0.04);

  b = @(T) 0.5*(r(T)-l(T));
  db = @(T) 0.5*(dr(T)-dl(T));

  c = @(T) 0.5*(r(T)+l(T));
  dc = @(T) 0.5*(dr(T)+dl(T));

  sq = @(T,x) sqrt((x-m(T)).^2 + q(T).^2);
  dsq = @(T,x)  ( -(x-m(T))*dm + q(T).*dq(T) )./sq(T,x);
  %
  % The total implied variance function w(T,x)
  %
  w = @(T,x) a(T) + b(T).*(x-m(T)) + c(T).*sq(T,x);

  dwdt = @(T,x) da(T) + db(T).*(x-m(T)) - b(T).*dm + ...
         dc(T).*sq(T,x) + c(T).*dsq(T,x);

  dwdx = @(T,x) b(T) + c(T).*(x-m(T))./sq(T,x);
  d2wdx2 = @(T,x) c(T)./sq(T,x) - c(T)*(x-m(T))*(x-m(T))./sq(T,x).^3;
  d2wdx2 = @(T,x) c(T).*q(T).^2./sq(T,x).^3;
  %
  % Dupires formula, but for SVI parametrization (the implied volatility is
  % approximated, not the option values). (Does not look exactly as in some
  % other descriptions.)
  %
  % The terms in the denominator first
  %
  n1 = @(T,x) (1-0.5*x.*dwdx(T,x)./w(T,x)).^2;
  n2 = @(T,x) -(0.5*dwdx(T,x)./w(T,x)).^2.*((0.5*w(T,x).*T + 1).^2 - 1);
  n3 = @(T,x) 0.5*T.*d2wdx2(T,x);

  w_local = @(T,x) (w(T,x) + T.*dwdt(T,x))./(n1(T,x)+n2(T,x)+n3(T,x));
  %
  % Take the input s values and convert to x
  %
  xx = log(ss./(s0*exp(ir*tt)));
  %
  % Evaluate the local volatility surface
  %
  sz = size(tt);
  pos = find(tt(:)==0);
  if (length(pos)>0)
    tt=tt(:);
    tt(pos) = eps;
    tt = reshape(tt,sz);
  end  
  ww = w_local(tt,xx);
  sigma = sqrt(ww);



