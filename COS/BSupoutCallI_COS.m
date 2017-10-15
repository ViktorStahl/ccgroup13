function [U] = BSupoutCallI_COS(S,K,T,r,sig,B)
% BENCHOP Problem 1: The Black-Scholes-Merton model for one underlying asset
% BSupoutCallI_COS computes the price for a barrier call up-and-out option (Standard parameters)
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

% Loop for Richardson extrapolation
MRc = 1;
for M = 2.^[4:7]

    % Parameters
    x = log(S./K);
    h = log(B/K);

    % Time step
    dt = T/M;
    
    % Interval [a,b]
    a = log(60/K); 
    b = log(140/K); 

    % Number of Fourier cosine coefficients
    N = 2^7;
    k = 0:N-1;
    omega = k'*pi/(b-a);
    temp = omega(2:end);
    temp = temp';
    
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
    
    % Characteristic function
    cf = exp((r-0.5*sig^2)*1i*omega*dt-0.5*sig^2*omega.^2*dt); 
    

    k=1:N-1;
    k=k';
    
    % Elements m_j
    mj=(exp(1i*(h-a)*temp')-1)./k;        
    mj_0=1i*pi*(h-a)/(b-a);        
    mj_N=(exp(1i*N*pi*(h-a)/(b-a))-1)/N;  
    mj_minus=-conj(mj); 
    mfactor1=exp(1i*N*pi*(h-a)/(b-a));        
    mfactor2=1;        
    mj_add=(mfactor1*exp(1i*(h-a)*temp')-mfactor2*1)./(k+N);
    
    % Vector m_s
    ms=[mj_0; mj_minus; 0; mj(end:-1:1)];
    % Vector m_c
    mc=[mj_add(end:-1:1); mj_N; mj(end:-1:1); mj_0];

    % Storage of vectors for the computation of the Fourier cosine
    % coefficients with the help of FFT algorithm
    fftms = fft(ms);
    fftmc = fft(mc);
    sgnvector = ones(2*N,1);
    sgnvector(2:2:end) = -1;
       
    if M > 1

        for m = M-1:-1:1

            % Vector u
            uj=cf'.*Uk;
            % Vector u_s
            us=[uj'; zeros(N,1)];           
            
            % Matrix-vector mulitplication M_s*u with the help of FFT algorithm
            fftu=fft(us);
            Msu=ifft(fftms.*fftu);
            Msu=Msu(1:N);
            % Matrix-vector mulitplication M_c*u with the help of FFT algorithm
            Mcu=ifft(fftmc.*sgnvector.*fftu);
            Mcu=Mcu(1:N);
            Mcu=Mcu(end:-1:1);

            % Fourier cosine coefficients option value Uk(t_m)
            Uk = (exp(-r*dt)/pi)*imag(Msu+Mcu)'; 
            Uk(1) = 0.5*Uk(1);
            
        end%m

    end%M>1
       
    % Fourier cosine coefficients density function 
    Recf = repmat(exp(-0.5*sig^2*dt*omega.^2),1,length(x)).*cos(omega*((r-0.5*sig^2)*dt+x-a));

    % Option value
    U = exp(-r*dt)*2/(b-a)*K*(Uk*Recf);

    U_Rich(MRc,:) = U;
    
    MRc = MRc+1;
end%M    

% 4-point Richardson extrapolation scheme (with k0=1/2, k1=1, k2=3/2)
U = 1/(5-3*sqrt(2))*(8*U_Rich(4:end,:)-(6*sqrt(2)+4)*U_Rich(3:end-1,:)+(3*sqrt(2)+2)*U_Rich(2:end-2,:)-U_Rich(1:end-3,:));

end