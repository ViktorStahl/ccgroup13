function [Delta]=BSeuCallDeltaI_COS(S,K,T,r,sig)
% BENCHOP Problem 1: The Black-Scholes-Merton model for one underlying asset
% BSeuCallDeltaI_COS computes the Delta for a European call option (Standard parameters)
%
% Input:    S       - Initial asset price   
%           K       - Strike price
%           T       - Terminal time  
%           r       - Risk-free interest rate
%           sig     - Volatility
%
% Output:   Delta   - Option's Delta
%
% This MATLAB code has been written for the BENCHOP project and is based on 
% the COS methodes developed by F. Fang, C.W. Oosterlee, and M.J. Ruijter
% Copyright 2015 by M.J. Ruijter

% Parameters
x = log(S./K);

% Cumulants
c1 = (r-0.5*sig^2)*T;
c2 = sig^2*T;
% Interval [a,b]
L = 8;
a = c1-L*sqrt(c2);
b = c1+L*sqrt(c2);

% Number of Fourier cosine coefficients
N = 20;
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

% Fourier cosine coefficients density function
Recf_Delta = repmat(exp(-0.5*sig^2*T*omega.^2).*-omega,1,length(x)).*sin(omega*((r-0.5*sig^2)*T+x-a));

% Option's Delta
Delta = exp(-r*T)*2/(b-a)*K*(Uk*Recf_Delta)./S;

end