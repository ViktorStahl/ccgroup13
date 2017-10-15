function [U] = MRTeuCall_COS(S,K,T,r,sig,lambda,gamma,delta)
% BENCHOP Problem 5: The Merton jump diffusion model for one underlying asset
% MRTeuCall_COS computes the price for a European call option
%
% Input:    S       - Initial asset price   
%           K       - Strike price
%           T       - Terminal time  
%           r       - Risk-free interest rate
%           sig     - Volatility
%           lambda  - Intensity rate Poisson process 
%           gamma   - Mean of log-jump size
%           delta   - Volatility of log-jump size
%
% Output:   U       - Option value
%
% This MATLAB code has been written for the BENCHOP project and is based on 
% the COS methodes developed by F. Fang, C.W. Oosterlee, and M.J. Ruijter
% Copyright 2015 by M.J. Ruijter

% Parameters
x = log(S./K);
xi = exp(gamma+0.5*delta^2)-1;

% Cumulants
c1 = (r-lambda*xi-0.5*sig^2+lambda*gamma)*T; 
c2 = (sig^2+lambda*gamma^2+delta^2*lambda)*T;
c4 = lambda*(gamma^4+6*delta^2*gamma^2+3*delta^4*lambda)*T;
% Interval [a,b]
L = 6;
a = c1-L*sqrt(c2+sqrt(c4));   
b = c1+L*sqrt(c2+sqrt(c4));

% Number of Fourier cosine coefficients
N = 70;
k = 0:N-1;
omega = k'*pi/(b-a);

% Fourier cosine coefficients payoff function
omegamina = omega*(-a);
omegabmina = omega*(b-a);
sin_omegamina = sin(omegamina);
sin_omegabmina = sin(omegabmina);
cos_omegamina = cos(omegamina);
cos_omegabmina = cos(omegabmina);
chi = ((cos_omegabmina+omega.*sin_omegabmina)*exp(b) ...
      -cos_omegamina -omega.*sin_omegamina)./(1+omega.^2);
psi = (sin_omegabmina-sin_omegamina)./omega;
psi(1) = b;
Uk = (chi-psi)';
Uk(1) = 0.5*Uk(1);

% Characteristic function
phiJ = exp(1i*gamma*omega-0.5*delta^2*omega.^2);
cf = exp((r-lambda*xi-0.5*sig^2)*1i*omega*T-0.5*sig^2*omega.^2*T).*exp(lambda*T*(phiJ-1));

% Fourier cosine coefficients density function
Recf = real(repmat(cf,1,length(x)).*exp(1i*omega*(x-a)));

% Option value 
U = exp(-r*T)*2/(b-a)*K*(Uk*Recf);

end