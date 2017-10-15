function [U]=BSeuCallUII_COS(S,K,T,r,sig)
% BENCHOP Problem 1: The Black-Scholes-Merton model for one underlying asset
% BSeuCallUII_COS computes the price for a European call option (Challenging parameters)
%
% Input:    S       - Initial asset price   
%           K       - Strike price
%           T       - Terminal time  
%           r       - Risk-free interest rate
%           sig     - Volatility
%
% Output:   U       - Option value
%
% This MATLAB code has been written for the BENCHOP project and is based on 
% the COS methodes developed by F. Fang, C.W. Oosterlee, and M.J. Ruijter
% Copyright 2015 by M.J. Ruijter

% Parameters
x = log(S./K);

% Interval [a,b]
a = log(60/K); 
b = log(170/K);

% Number of Fourier cosine coefficients
N = 234;
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
Recf = repmat(exp(-0.5*sig^2*T*omega.^2),1,length(x)).*cos(omega*((r-0.5*sig^2)*T+x-a));

% Option value
U = exp(-r*T)*2/(b-a)*K*(Uk*Recf);

end