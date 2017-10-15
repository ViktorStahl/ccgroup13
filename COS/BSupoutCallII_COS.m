function [U] = BSupoutCallII_COS(S,K,T,r,sig,B)
% BENCHOP Problem 1: The Black-Scholes-Merton model for one underlying asset
% BSupoutCallII_COS computes the price for a barrier call up-and-out option (Challenging parameters)

%
% Input:    S       - Initial asset price   
%           K       - Strike price
%           T       - Terminal time  
%           r       - Risk-free interest rate
%           sig     - Volatility
%           B       - Barrier
%
% Output:   U       - Option value
%
% This MATLAB code has been written for the BENCHOP project and is based on 
% the COS methodes developed by F. Fang, C.W. Oosterlee, and M.J. Ruijter
% Copyright 2015 by M.J. Ruijter

% Parameters
x = log(S./K);
h = log(B/K);

% Time step
dt = T;

% Interval [a,b]
a = log(60/K); 
b = log(128/K); 

% Number of Fourier cosine coefficients
N = 187;
k = 0:N-1;
omega = k'*pi/(b-a);

% Fourier cosine coefficients payoff function
omegamina = omega*(-a);
omegahmina = omega*(h-a);
sin_omegamina = sin(omegamina);
sin_omegahmina = sin(omegahmina);
cos_omegamina = cos(omegamina);
cos_omegahmina = cos(omegahmina);
chi = ((cos_omegahmina+omega.*sin_omegahmina)*exp(h)...
      -cos_omegamina-omega.*sin_omegamina)./(1+omega.^2);
psi = (sin_omegahmina-sin_omegamina)./omega;
psi(1) = h;  
Uk = (chi-psi)';
Uk(1) = 0.5*Uk(1);

% Fourier cosine coefficients density function
Recf = repmat(exp(-0.5*sig^2*dt*omega.^2),1,length(x)).*cos(omega*((r-0.5*sig^2)*dt+x-a));

% Option value
U = exp(-r*dt)*2/(b-a)*K*(Uk*Recf);

end