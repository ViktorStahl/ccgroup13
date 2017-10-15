function [U] = BSeuLocVolI_COS(S,K,T,r,sig)
% BENCHOP Problem 3: The Black-Scholes-Merton model with local volatility
% BSeuLocVolI_COS computes the price for a European call option
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

% Number of timesteps
M = 17;
% Timestep
dt = T/M;

% Drift and volatility functions and their derivatives
fmu = @(S)r*S;   
fmu_S = @(S)r;
fsigma = @(t,S)(0.15+0.15.*(0.5+(2.*t)).*(S./0.1e3-1.2).^2./(S.^2./0.1e5+0.144e1)).*S;
fsigma_S = @(t,S)((0.225e-8+0.3e-8.*t).*S.^4+(0.648e-4+0.864e-4.*t).*S.^2+(-0.5184e-2-0.20736e-1.*t).*S+0.46656+0.62208.*t)./(0.144e1+S.^2./0.1e5).^2;
fsigma_SS = @(t,S)(0.5+(2.*t)).*(0.31104e-5.*S.^2-0.1492992e-1)./(0.144e1+S.^2./0.1e5).^3;
fsigma_t = @(t,S)0.3.*(-0.12e1+S./0.1e3).^2.*S./(0.144e1+S.^2./0.1e5);

% Interval [a,b]
a = 50; b = 360;

% Number of Fourier cosine coefficients
N = 2^5;
k = 0:N-1;
omega = k'*pi/(b-a);

% Fourier cosine coefficients payoff function  
omegaKmina = omega*(K-a);
omegabmina = omega*(b-a);
sin_omegaKmina = sin(omegaKmina);
sin_omegabmina = sin(omegabmina);
cos_omegaKmina = cos(omegaKmina);
cos_omegabmina = cos(omegabmina);
ksi = (cos_omegabmina+omega.*(-K*sin_omegaKmina+1/(b-a)*sin_omegabmina*b) ...
      -cos_omegaKmina)./(omega.^2);
ksi(1) = b^2/2-K^2/2;
psi = (sin_omegabmina-sin_omegaKmina)./omega;
psi(1) = b-K;
Uk = 2/(b-a)*(ksi-K*psi)';
Uk(1) = 0.5*Uk(1);

% Grid parameters and grid Sgrid
dS = (b-a)/N;
Sgrid = (a+0.5*dS):dS:(b-0.5*dS);
expiomegaSmina = exp(1i*omega*(Sgrid-a));

% Drift and its derivatives on Sgrid
mu_Sgrid = fmu(Sgrid);
mu_S_Sgrid = fmu_S(Sgrid);
% Drift and its derivatives on S
mu_S = fmu(S);
mu_S_S = fmu_S(S);

for mm = M-1:-1:1

    % Time
    tm = mm*dt;
    
    % Volatility and its derivatives on Sgrid at time tm
    sigma = fsigma(tm,Sgrid);
    sigma_S = fsigma_S(tm,Sgrid);
    sigma_SS = fsigma_SS(tm,Sgrid);
    sigma_t = fsigma_t(tm,Sgrid);
    
    % Characteristic function at time tm
    m = mu_Sgrid-0.5*sigma.*sigma_S+0.5*mu_Sgrid.*mu_S_Sgrid*dt;
    s = sigma + 0.5*(sigma_t+mu_S_Sgrid.*sigma+mu_Sgrid.*sigma_S+0.5*sigma_SS.*sigma.^2)*dt;
    kappa = 0.5*sigma.*sigma_S;
    onemin2omegakappadt=(1-2*1i*omega*kappa*dt);
    cf = exp(1i*omega*m*dt).*exp(-0.5*(omega.^2)*(s.^2)*dt./onemin2omegakappadt).*onemin2omegakappadt.^(-1/2); 
    
    % Fourier cosine coefficients density function
    Recf = real(cf.*expiomegaSmina);
    
    % Option value   
    U = exp(-r*dt)*(Uk*Recf);

    % Fourier cosine coefficients option value Uk(t_m) 
    % with the help of DCT
    Uk = dct(2/(b-a)*U*dS);
    Uk(2:end) = 1/sqrt(2/N)*Uk(2:end);
    Uk(1) = 1/sqrt(1/N)*Uk(1);
    Uk = Uk(1:N); 
    Uk(1) = 0.5*Uk(1);
    
end%mm

% Volatility and its derivatives on S at time 0
sigma = fsigma(0,S);
sigma_S = fsigma_S(0,S);
sigma_SS = fsigma_SS(0,S);
sigma_t = fsigma_t(0,S);

% Characteristic function at time 0
m = mu_S-0.5*sigma.*sigma_S+0.5*mu_S.*mu_S_S*dt;
s = sigma + 0.5*(sigma_t+mu_S_S.*sigma+mu_S.*sigma_S+0.5*sigma_SS.*sigma.^2)*dt;
kappa = 0.5*sigma.*sigma_S;
onemin2omegakappadt=(1-2*1i*omega*kappa*dt);
cf = exp(1i*omega*m*dt).*exp(-0.5*(omega.^2)*(s.^2)*dt./onemin2omegakappadt).*onemin2omegakappadt.^(-1/2); 

% Fourier cosine coefficients density function
Recf = real(cf.*exp(1i*omega*(S-a)));

% Option value   
U = exp(-r*dt)*(Uk*Recf);

end