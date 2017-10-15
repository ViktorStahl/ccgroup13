function [U]=BSamCallDD_COS(S,K,T,r,sig,D,alpha)
% BENCHOP Problem 2: The Black-Scholes-Merton model with discrete dividends
% BSamCallDD_COS computes the price for an American call option
%
% Input:    S       - Initial asset price   
%           K       - Strike price
%           T       - Terminal time  
%           r       - Risk-free interest rate
%           sig     - Volatility
%           D       - At time tau = alpha*T a dividend D*S is paid
%           alpha   - At time tau = alpha*T a dividend D*S is paid
%
% Output:   U       - Option value
%
% This MATLAB code has been written for the BENCHOP project and is based on 
% the COS methodes developed by F. Fang, C.W. Oosterlee, and M.J. Ruijter
% Copyright 2015 by M.J. Ruijter

% Parameters
x = log(S./K);
tau = alpha*T;

% Cumulants
c1 = (r-0.5*sig^2)*T;
c2 = sig^2*T;
% Interval [a,b]
L = 8;
a = c1-L*sqrt(c2);
b = c1+L*sqrt(c2);

% Number of Fourier cosine coefficients
N = 137;
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
Uk = 2/(b-a)*K*(chi-psi)';
Uk(1) = 0.5*Uk(1);

% Grid ymin and yplus
dy = (b-a)/N;
ymin = (a+0.5*dy):dy:(b-0.5*dy);
yplus = ymin+log(1-D);

% Fourier cosine coefficients density function (time-interval [tau,T])
Recf = repmat(exp(-0.5*sig^2*(T-tau)*omega.^2),1,length(yplus)).*cos(omega*((r-0.5*sig^2)*(T-tau)+yplus-a));

% Continuation value and payoff (time tau-)
c = exp(-r*(T-tau))*Uk*Recf;
g = K*(exp(ymin)-1);
% Option value (time tau-)
U = max(c,g);
% Fourier cosine coefficients option value Uk(tau-) (time tau-)
Uk = dct(2/(b-a)*U*dy);
Uk(2:end) = 1/sqrt(2/N)*Uk(2:end);
Uk(1) = 1/sqrt(1/N)*Uk(1);
Uk = Uk(1:N);     
Uk(1) = 0.5*Uk(1);

% Fourier cosine coefficients density function (time-interval [0,tau])
Recf = repmat(exp(-0.5*sig^2*tau*omega.^2),1,length(x)).*cos(omega*((r-0.5*sig^2)*tau+x-a));

% Option value (time 0)
U = exp(-r*tau)*Uk*Recf;

end