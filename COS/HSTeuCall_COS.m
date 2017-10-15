function [U] = HSTeuCall_COS(S,K,T,r,V,kap,th,sig,rho)
% BENCHOP Problem 4: The Heston model for one underlying asset
% HSTeuCall_COS computes the price for a European call option
%
% Input:    S       - Initial asset price   
%           K       - Strike price
%           T       - Terminal time  
%           r       - Risk-free interest rate
%           V       - Initial instantaneous variance
%           kap     - Speed of mean reversion
%           th      - Long run variance
%           sig     - Volatility of variance process
%           rho     - Correlation coefficient
%
% Output:   U       - Option value
%
% This MATLAB code has been written for the BENCHOP project and is based on 
% the COS methodes developed by F. Fang, C.W. Oosterlee, and M.J. Ruijter
% Copyright 2015 by M.J. Ruijter

% Parameters
x = log(S./K);

% Cumulants
c1 = r*T+(1-exp(-kap*T))*(th-V)/(2*kap)-0.5*th*T;
c2alt = th*(1+sig)*T;
% Interval [a,b]
L = 8;
a = c1-L*abs(sqrt(c2alt));
b = c1+L*abs(sqrt(c2alt));

% Number of Fourier cosine coefficients
N = 28;
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
D = sqrt((kap-1i*rho*sig*omega).^2+(omega.^2+1i*omega)*sig^2);
G = (kap-1i*rho*sig*omega-D)./(kap-1i*rho*sig*omega+D);
cf = exp(1i*omega*r*T+V/sig^2*(1-exp(-D*T))./(1-G.*exp(-D*T)) ...
     .*(kap-1i*rho*sig*omega-D))...
     .*exp(kap*th/sig^2*(T*(kap-1i*rho*sig*omega-D)-2*log((1-G.*exp(-D*T))./(1-G))));

% Fourier cosine coefficients density function
Recf = real(repmat(cf,1,length(x)).*exp(1i*omega*(x-a)));

% Option value 
U = exp(-r*T)*2/(b-a)*K*(Uk*Recf);

end