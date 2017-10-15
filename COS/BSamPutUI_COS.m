function [U] = BSamPutUI_COS(S,K,T,r,sig)
% BENCHOP Problem 1: The Black-Scholes-Merton model for one underlying asset
% BSamPutUI_COS computes the price for an American put option (Standard parameters)
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

% Loop for Richardson extrapolation
MRc = 1;
for M = 2.^[2:5]
   
    % Parameters
    x = log(S./K);

    % Time step
    dt = T/M;
    
    % Interval [a,b]
    a = log(50/K); 
    b = log(160/K); 

    % Number of Fourier cosine coefficients
    N = 2^6;
    k = 0:N-1;
    omega = k'*pi/(b-a);
    temp = omega(2:end);
    temp = temp';

    % Fourier cosine coefficients payoff function
    omegamina = omega*(-a);
    sin_omegamina = sin(omegamina);
    cos_omegamina = cos(omegamina);
    chi = (cos_omegamina-exp(a)+omega.*sin_omegamina)./(1+omega.^2);
    psi = sin_omegamina./omega;
    psi(1) = -a;    
    Uk = (-chi+psi)';
    Uk(1) = 0.5*Uk(1);
    
    % Characteristic function
    cf = exp((r-0.5*sig^2)*1i*omega*dt-0.5*sig^2*omega.^2*dt); 
    partcf = exp(-0.5*sig^2*dt*omega.^2);

    % xs is the early-exercise point where c = g,
    xs = 0; % initial value

    if M > 1

        for m = M-1:-1:1

            % Newton-Raphson iterations to find the early-exercise point
            % where f = c-g = 0
            for NR = 1:5

                % Fourier cosine coefficients density function
                Recf = partcf.*cos(omega*((r-0.5*sig^2)*dt+xs-a));
                Recf_x = partcf.*-omega.*sin(omega*((r-0.5*sig^2)*dt+xs-a));                
                
                % Continuation and payoff value and derivatives
                c = exp(-r*dt)*2/(b-a)*K*(Uk*Recf);
                c_x = exp(-r*dt)*2/(b-a)*K*(Uk*Recf_x);
                g = 0;
                g_x = 0;

                if xs <= 0   
                    g = K*(1-exp(xs));  % x = log(S/K),S = Kexp(x)
                    g_x = -K*exp(xs);
                end

                f = c-g;
                f_x = c_x-g_x;

                % Next approximation to the root
                xs = xs-f/f_x;

            end%NR
            if xs < a
                xs = a;
            elseif xs > b
                xs = b;
            end

            % Fourier cosne coefficients payoff function Gk(a,x2)            
            omegaxsmina = omega*(xs-a);
            sin_omegaxsmina = sin(omegaxsmina);
            cos_omegaxsmina = cos(omegaxsmina);
            chi = (cos_omegaxsmina*exp(xs)-exp(a)+omega.*sin_omegaxsmina*exp(xs))./(1+omega.^2);
            psi = sin_omegaxsmina./omega;
            psi(1) = xs-a;
            G = (-chi+psi)';
            % Fourier cosine coefficients continuation value Ck(t_m,xs,b)
            C = Cvalue(xs,b,N,a,b,temp,cf,Uk,dt,r);
            % Fourier cosine coefficients option value Uk(t_m)
            Uk = C+G;
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

% 4-point Richardson extrapolation scheme (with k0=1, k1=2, k2=3)
U = 1/21*(64*U_Rich(4:end,:)-56*U_Rich(3:end-1,:)+14*U_Rich(2:end-2,:)-U_Rich(1:end-3,:));

end


function C=Cvalue(x1, x2, N, a, b, temp, cf, U, dt, r)
%
% Input:    x1      - Lower bound integral
%           x2      - Upper bound integral
%           N       - Number of terms series expansion
%           a       - Left bound computational domain
%           b       - Right bound computational domain   
%           temp    - Vector with values [1:N-1]*pi/(b-a);
%           cf      - Characteristic function
%           U       - Fourier cosine coefficients Uk(t_{m+1})
%           dt      - Time-step
%           r       - Risk-free interest rate
%
% Output:   C       - Fourier cosine coefficients Ck(t_m,x1,x2)

k=1:N-1;
k=k';
         
% Elements m_j
mj=(exp(1i*(x2-a)*temp')-exp(1i*(x1-a)*temp'))./k;        
mj_0=1i*pi*(x2-x1)/(b-a);        
mj_N=(exp(1i*N*pi*(x2-a)/(b-a))-exp(1i*N*pi*(x1-a)/(b-a)))/N;  
mj_minus=-conj(mj); 
mfactor1=exp(1i*N*pi*(x2-a)/(b-a));        
mfactor2=exp(1i*N*pi*(x1-a)/(b-a));        
mj_add=(mfactor1*exp(1i*(x2-a)*temp')-mfactor2*exp(1i*(x1-a)*temp'))./(k+N);
  
% Vector m_s
ms=[mj_0; mj_minus; 0; mj(end:-1:1)];
% Vector m_c
mc=[mj_add(end:-1:1); mj_N; mj(end:-1:1); mj_0];

% Vector u
uj=cf'.*U;
% Vector u_s
us=[uj'; zeros(N,1)];
        
% Matrix-vector mulitplication M_s*u with the help of FFT algorithm
fftu=fft(us);
Msu=ifft(fft(ms).*fftu);
Msu=Msu(1:N);
% Matrix-vector mulitplication M_c*u with the help of FFT algorithm
sgnvector=ones(2*N,1);
sgnvector(2:2:end)=-1;
Mcu=ifft(fft(mc).*sgnvector.*fftu);
Mcu=Mcu(1:N);
Mcu=Mcu(end:-1:1);

% Fourier cosine coefficients Ck(t_m,x1,x2)
C=(exp(-r*dt)/pi)*imag(Msu+Mcu); 
C=C';

end