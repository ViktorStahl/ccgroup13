function [S V_out payoff] = FDMUniformGrid( is_call, is_American_style, ...
    T, X, barrier_type, H, ...
    r, q, volatility_model, volatility_parameter, ...
    div_datetimes, div_amounts, div_percentages, ... 
    n_price_steps_per_strike, Smax_multiplier, stick_strike_to_node, ...
    n_time_steps_per_year, min_time_steps_per_interval, n_Rannacher_sub_steps)
%
% FDMUniformGrid solves Black-Scholes PDE for financial option price with multiple discrete dividends 
% by Crank-Nicolson scheme with Rannacher time stepping. The grid in underlying price S 
% is uniform between lower and upper boundaries. Time steps are uniform over each interval 
% of continuity (between maturity T, successive dividend dates and current date t=0).
%
% Author: Yuri Shpolyanskiy, Tbricks AB, www.tbricks.com, University ITMO, en.ifmo.ru/en/
% Date: November 2014
%
% References:
% 1. P.Willmott, On quantitative finance, Wiley, 2nd ed. (2006), 1500 p, Chapter 78.
% 2. R.Rannacher, Finite element solution of diffusion problems with irregular data, 
%    Numerische Mathematik 43, N2 (1984), pp. 309-327. 
% 3. S.Ikonen, J.Toivanen, Pricing American options using LU decomposition, 
%    Appl. Math. Sci. 1, N51 (2007), pp. 2529-2551.
%
% function [S V_out payoff] = FDMUniformGrid( is_call, is_American_style, ...
%     T, X, barrier_type, H, ...
%     r, q, volatility_model, sigma, ...
%     div_datetimes, div_amounts, div_percentages, ... 
%     n_price_steps_per_strike, Smax_multiplier, stick_strike_to_node, ...
%     n_time_steps_per_year, min_time_steps_per_interval, n_Rannacher_sub_steps)
%
% Input arguments:
% is_call: 1 for call option, 0 for put option
% is_American_style: 1 for American options, 0 for European options
% T: time to maturity in years
% X: strike price
% barrier_type: option type, 0 for vanilla option, 1 for up-and-out barrier
% H: barrier level for barrier options
% r: annualized interest rate
% q: annualized dividend yield
% volatility_model: volatility model, 1 for constant volatility
% volatility_parameter: annualized constant volatility if volatility_model==1
%                       cell array if volatility_model>1:
%                       {handle_to_local_vol_function, S0} 
% div_datetimes: vector of dividend dates
% div_amounts: vector of dividend amounts for absolute discrete dividends
% div_percentages: vector of dividend ratios for proportional discrete dividends
% n_price_steps_per_strike: number of price steps per strike
% Smax_multiplier: multiplier defining upper level Smax of underlying price for vanilla options as
%                  Smax = X * Smax_multiplier
% stick_strike_to_node: 1 to guarantee that strike prices coincides with one of price nodes, 
%                       0 to have strike price between price nodes
% n_time_steps_per_year: requested number of time steps per year
% min_time_steps_per_interval: minimum number of time steps per interval of continuity
% n_Rannacher_sub_steps: number of sub-steps of Backward Euler scheme in the beginning of
%                        each interval of continuity. First time step dt in the interval is replaced by 
%                        n_Rannacher_sub_steps sub-steps of duration dt/n_Rannacher_sub_steps.
% 
% Output arguments:
% S: grid of underlying price levels
% V_out: grid of option prices corresponding to S
% payoff: grid of payoff values corresponding to S



if X < 1e-10
    error('Strike price should be positive');
end

if T < 0
    error('Time to maturity cannot be negative');
end

if is_American_style && barrier_type > 0
    error('American-style barrier options are not supported')
end

sigma = 0;
if volatility_model == 1
    sigma = volatility_parameter;
    if sigma < 1e-4
        error('Constant volatility cannot be lower than 1e-4')
    end
else
    locvolfunc = cell2mat( volatility_parameter(1) );
    S0 = cell2mat( volatility_parameter(2) );
end

if n_price_steps_per_strike < 20
    error('Number of price steps should not be less than 20');
end    

if n_time_steps_per_year < 20
    error('Number of time steps per year cannot be less than 20');
end

if n_Rannacher_sub_steps < 0
    error('Number of Rannacher sub-steps cannot be negative');
end

switch barrier_type
    case 0 % Vanilla option
        if Smax_multiplier < 1.1999999999 || Smax_multiplier > 12.0000000001    
            error('Smax_multiplier should be between 1.2 and 12.0');        
        end    
        Smin = 0;    

        n_price_steps = ceil( Smax_multiplier * n_price_steps_per_strike );
        if stick_strike_to_node
            dS_approx = X / n_price_steps_per_strike;
        else
            dS_approx = X / (n_price_steps_per_strike - 0.5);
        end
        Smax = Smin + n_price_steps * dS_approx;
        
    case 1 %Barrier up and out option
        if ~is_call
            error('Up-and-out put options are not supported yet')
        end
        if (H <= X)
            error('Barrier level should be higher than the strike price for up-and-out call option')
        end
        Smin = 0;
        Smax = H;
        if stick_strike_to_node
            n_price_steps = round(H / X * n_price_steps_per_strike);
        else
            dS_approx = X / (n_price_steps_per_strike - 0.5);
            n_price_steps = round(H / dS_approx);
        end
        
    otherwise
        error('Unsupported type of barrier option');
end

% Minimum and maximum indices of unknowns for differnt boundary conditions
bc_Smin_Dirichlet = 1;

bc_Smax_Dirichlet = 1;
if is_call && barrier_type == 0
    bc_Smax_Dirichlet = 0;
end

n_price_levels = n_price_steps+1; 
S = linspace(Smin, Smax, n_price_levels); 
dS = S(2) - S(1);

payoff(n_price_levels) = 0.0; %preallocation
if is_call
    i_X = ceil((X - Smin) / dS) + 1;
    payoff(i_X:n_price_levels) = S(i_X:n_price_levels) - X;
else    
    i_X = floor((X - Smin) / dS);
    payoff(1:i_X) = X - S(1:i_X);
end

if bc_Smax_Dirichlet
    payoff(n_price_levels) = 0.0;
end

nD = length(div_datetimes);
if nD > 0
    nD_absolute = length(div_amounts);
    if nD_absolute > 0 && nD_absolute ~= nD
        error('Wrong size of vector with absolute dividends. It should be either zero or equal to the number of dividend datetimes.')
    end
    nD_percentage = length(div_percentages);
    if nD_percentage > 0 && nD_percentage ~= nD
        error('Wrong size of vector with percentage dividends. It should be either empty or equal to the number of dividend datetimes.')
    end
end

for iD = 1 : nD
    if div_datetimes(iD) < 0
        error('Dividend datetimes cannot be negative')
    end
    if iD>1
        if div_datetimes(iD) <= div_datetimes(iD-1)
            error('Dividend datetimes should be sorted in ascending order')
        end
    end
end

V_out = payoff;
t_in = T;

iD = nD;

if iD > 0
    tD = div_datetimes(iD);
    if tD > T + 1e-10
        error('Last dividend datetime exceeds maturity')
    end
    if tD > T - 1e-8
        V_in = V_out;
        ApplyDividendMinPercentageAbsolute(V_in, iD);
        t_in = div_datetimes(iD);
        iD = iD - 1;
    end
end

if T < 1e-8
    return;
end

scheme_CN = 1; %constant value for nested functions
scheme_BE = 2; %constant value for nested functions

scheme = 0; %variable to keep current FDM scheme for nested functions

A(n_price_levels) = 0.0;     %lower diagonal of tridiagonal system
B(n_price_levels) = 0.0;     %main diagonal of tridiagonal system
C(n_price_levels) = 0.0;     %upper diagonal of tridiagonal system
RHS(n_price_levels) = 0.0;   %right-hand-side of tridiagonal system 
u(n_price_levels) = 0.0;     %vector used in LU decomposition
l(n_price_levels) = 0.0;     %vector used in LU decomposition
inv_d(n_price_levels) = 0.0; %vector used in LU decomposition
w(n_price_levels) = 0.0;     %vector used in LUSolve

i_exercise = -1; %index of exercise boundary

if bc_Smin_Dirichlet
    i_first = 2; %Dirichlet boundary condition at Smin
else
    error('Unsupported BC type at Smin');
end

if bc_Smax_Dirichlet
    i_last = n_price_levels - 1; %Dirichlet boundary condition at Smax
else
    i_last = n_price_levels; %Zero second derivative at Smax
end

need_LU_decompose = 1;

t = T; %current time to be used in nested functions
dt = 0; %variable to keep time step will be used in nested functions

while iD > 0
    t_out = div_datetimes(iD);
    V_in = V_out;
    RunInterval(t_in, t_out, V_in);
  
    V_in = V_out;
    ApplyDividendMinPercentageAbsolute(V_in, iD);

    iD = iD - 1;
    t_in = t_out;
end

t_out = 0;
V_in = V_out;
RunInterval(t_in, t_out, V_in);

%NESTED FUNCTONS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = RunInterval(t_in, t_out, V_in)

    t_interval = t_in - t_out;

    V_out = V_in;
    if abs(t_interval) < 1e-8
        return;
    end

    if t_out > t_in
        error('t_out should be less than t_in in due to backward time direction');
    end
    
    t = t_in;

    dt_approx = 1.0 / n_time_steps_per_year;

    Nt = ceil( t_interval / dt_approx );
    if (Nt < min_time_steps_per_interval)
        Nt = min_time_steps_per_interval;
    end 

    dt_CN = t_interval / Nt;

    if n_Rannacher_sub_steps > 0
        SetScheme( scheme_BE, dt_CN / n_Rannacher_sub_steps );
        for s = 1 : n_Rannacher_sub_steps
            V_t = V_out;
            MakeBEStep(V_t); %MakeBEStep puts the result to V_out
        end
        n_lower = 2;
    else
        n_lower = 1;
    end
    
    SetScheme( scheme_CN, dt_CN );
    for n = n_lower : Nt
        V_t = V_out;
        MakeCNStep(V_t); %MakeCNStep puts the result to V_out
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = SetScheme( scheme_in, dt_in )
    scheme = scheme_in;
    dt = dt_in;
    %need_LU_decompose = 1; %it is set in CalculateABCVectors
    if scheme == scheme_CN
        CalculateABCVectorsCN();
    elseif scheme == scheme_BE
        CalculateABCVectorsBE();       
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = SolveTridiagonal( V_t )
    if is_American_style
        if is_call
            if need_LU_decompose
                LUDecompose();
                need_LU_decompose = 0;
            end
            
            LUSolveAmericanCall( V_t );
        else
            if need_LU_decompose
                ULDecompose();
                need_LU_decompose = 0;
            end
            
            ULSolveAmericanPut( V_t );
        end
    else    
        if need_LU_decompose
            LUDecompose();
            need_LU_decompose = 0;
        end
    
        LUSolve( V_t );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = LUDecompose()
    bc_diag_first=0;
    bc_superdiag_first=0;
    
    bc_diag_last=0;
    bc_subdiag_last=0;
    
    if ~bc_Smax_Dirichlet
        bc_diag_last      = (-2.0) * C(i_last);
        bc_subdiag_last   =          C(i_last);
    end
    
    i = i_first;
    inv_d(i) = 1.0 / (1.0 - B(i) + bc_diag_first);
    
    i = i+1;
    superdiag_prev_i = bc_superdiag_first - C(i-1); 
    u(i-1)   = superdiag_prev_i;
    l(i)     = - A(i) * inv_d(i-1);  
    inv_d(i) = 1.0 / ( 1.0 - B(i) - l(i) * superdiag_prev_i );
    
    for i = i_first + 2 : i_last - 1
        superdiag_prev_i = -C(i-1); 
        u(i-1)   = superdiag_prev_i;
        l(i)     = - A(i) * inv_d(i-1);  
        inv_d(i) = 1.0 / ( 1.0 - B(i) - l(i) * superdiag_prev_i );
    end
    
    i = i_last;
    superdiag_prev_i = -C(i-1); 
    u(i-1)   = superdiag_prev_i;
    l(i)     = ( bc_subdiag_last - A(i) ) * inv_d(i-1);  
    inv_d(i) = 1.0 / ( bc_diag_last + 1.0 - B(i) - l(i) * superdiag_prev_i );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = ULDecompose()
    bc_diag_first=0;
    bc_superdiag_first=0;
    
    bc_diag_last=0;
    bc_subdiag_last=0;

    if ~bc_Smax_Dirichlet
        bc_diag_first      = (-2.0) * A(i_first);
        bc_superdiag_first =          A(i_first);
    end
    
    i = i_last;
    
    inv_d(i) = 1.0 / ( 1.0 - B(i) + bc_diag_last );
    l(i) = - A(i) + bc_subdiag_last;
    
    for i = i_last - 1 : -1 : i_first + 1
        l(i) = - A(i);
        u(i) = - C(i) * inv_d(i+1);
        inv_d(i) = 1.0 / ( 1.0 - B(i) - l(i+1)*u(i) );
    end
    
    i = i_first;
    u(i) = ( -C(i) + bc_superdiag_first ) * inv_d(i+1);
    inv_d(i) = 1.0 / (1.0 - B(i) + bc_diag_first - l(i+1) * u(i) );
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LUSolve ( V_t )
    w = RHS;
    V_out = V_t;
    
    for i = i_first + 1 : i_last
        w(i) = w(i) - l(i) * w(i-1);
    end
    
    V_out(i_last) = w(i_last) * inv_d(i_last);
    for i = i_last - 1 : -1 : i_first
        V_out(i) = ( w(i) - u(i) * V_out(i+1) ) * inv_d(i);
    end
    
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LUSolveAmericanCall ( V_t )
    w = RHS;
    V_out = V_t;
    
    for i = i_first + 1 : i_last
        w(i) = w(i) - l(i) * w(i-1);
    end
    
    i = i_last;
    V_out(i_last) = w(i_last) * inv_d(i_last);
    
    if V_out(i_last) < payoff(i_last)
        V_out(i_last) = payoff(i_last);
        i_exercise = i_last;
        for i = i_last - 1 : -1 : i_first
            V_out(i) = ( w(i) - u(i) * V_out(i+1) ) * inv_d(i);
            if V_out(i) < payoff(i)
                V_out(i) = payoff(i);
            else
                i_exercise = i + 1;
                break;
            end
        end
    end
    
    i = i - 1;
    while i >= i_first
        V_out(i) = ( w(i) - u(i) * V_out(i+1) ) * inv_d(i);
        i = i - 1;
    end
    
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ULSolveAmericanPut ( V_t )
    w = RHS;
    V_out = V_t;
    
    for i = i_last - 1 : -1 : i_first
        w(i) = w(i) - u(i) * w(i+1);
    end
    
    i = i_first;
    V_out(i_first) = w(i_first) * inv_d(i_first);
    if V_out(i_first) < payoff(i_first)
        V_out(i_first) = payoff(i_first);
        i_exercise = i_first;
        for i = i_first + 1 : i_last
            V_out(i) = ( w(i) - l(i) * V_out(i-1) ) * inv_d(i);
            if V_out(i) < payoff(i)
                V_out(i) = payoff(i);
            else
                i_exercise = i - 1;
                break;
            end
        end
    end
    
    i = i + 1;
    while i <= i_last
        V_out(i) = ( w(i) - l(i) * V_out(i-1) ) * inv_d(i);
        i = i + 1;
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = CalculateABCVectorsCN()
%Calculate A, B, C vectors for CN scheme

    % Internal notation from Ref.1
    % nu1 = dt/dS^2; nu2 = dt/dS; 
    % ai = 0.5*(volatility*S(i))^2; bi = r*S(i);

    tmp_025_nu2 = 0.25 * dt / dS; 
    
    if volatility_model == 1
        tmp_05_nu1_ai  = tmp_025_nu2 / dS * sigma^2 * S.^2;
    elseif volatility_model == 2
        vol = locvolfunc(S, t); 
        tmp_05_nu1_ai  = tmp_025_nu2 / dS * ( S .* vol) .^ 2;
    elseif volatility_model == 3
        vol = real(locvolfunc(r, S0, S, t));
        tmp_05_nu1_ai  = tmp_025_nu2 / dS * ( S .* vol) .^ 2;
    else
        error(['Volatility model ' sprintf('%d',volatility_model) ' is not supported'] );
    end
    
    tmp_025_nu2_bi = tmp_025_nu2 * (r - q) * S;
    
    A =       tmp_05_nu1_ai  -  tmp_025_nu2_bi;
    B =(-2.0)*tmp_05_nu1_ai  -  0.5 * dt * r;
    C =       tmp_05_nu1_ai  +  tmp_025_nu2_bi;
        
    need_LU_decompose = 1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = CalculateABCVectorsBE()
%Calculate A, B, C vectors for BE scheme

    % Internal notation from Ref.1
    % nu1 = dt/dS^2; nu2 = dt/dS; 
    % ai = 0.5*(volatility*S(i))^2; bi = r*S(i);
    
    tmp_05_dt_div_dS = 0.5 * dt / dS; 
    
    if volatility_model == 1
        tmp_nu1_ai  = tmp_05_dt_div_dS / dS * sigma^2 * S.^2;
    elseif volatility_model == 2
        vol = locvolfunc(S, t); %LocalVolatility1(S, t);
        tmp_nu1_ai  = tmp_05_dt_div_dS / dS * ( S .* vol ) .^ 2;
    elseif volatility_model == 3
        vol = real(locvolfunc(r, S0, S, t));
        tmp_nu1_ai  = tmp_05_dt_div_dS / dS * ( S .* vol ) .^ 2;
    else
        error(['Volatility model ' sprintf('%d',volatility_model) ' is not supported'] );
    end
    
    tmp_05_nu2_bi = tmp_05_dt_div_dS * (r - q) * S;
    
    A =       tmp_nu1_ai  -  tmp_05_nu2_bi;
    B =(-2.0)*tmp_nu1_ai  -  dt * r   ;
    C =       tmp_nu1_ai  +  tmp_05_nu2_bi;
        
    need_LU_decompose = 1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = MakeCNStep( V_t )
    t_next = t - dt;

    V_out_Smin = 0;
    V_t_Smax = 0;
    if is_call
        %V_t_Smin = 0;
        %V_out_Smin = 0;
        
        if ~bc_Smax_Dirichlet
            %BC: Zero second derivative of V at Smax
            V_t_Smax = 2 * V_t(n_price_levels) - V_t(n_price_levels - 1);
        end
        
        %V_out_Smax = 0;
    else 
        %V_t_Smin = V_t(1);
        if is_American_style
            V_out_Smin = X;
        else
            V_out_Smin = X * exp( -r*(T-t_next));
        end
        
        %V_out_Smax = 0;
    end
    
    %for i = i_first : n_price_levels - 1
    for i = 2 : n_price_levels - 1
        RHS(i) = A(i) * V_t(i-1) + (B(i) + 1) * V_t(i) + C(i) * V_t(i+1);
    end
    
    if ~bc_Smax_Dirichlet %BC: Zero second derivative of V at Smax
        i = n_price_levels;
        RHS(i) = A(i) * V_t(i-1) + (B(i) + 1) * V_t(i) + C(i) * V_t_Smax;
    end
    
    t = t_next;
    
    if volatility_model > 1
        CalculateABCVectorsCN();
    end
    
    RHS(i_first) = RHS(i_first) + A(i_first) * V_out_Smin;
    %RHS(i_last ) = RHS(i_last ) + C(i_last)  * V_out_Smax;
    
    SolveTridiagonal( V_t );
    
    V_out(1) = V_out_Smin;
    
    if bc_Smax_Dirichlet
        V_out(n_price_levels) = 0;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = MakeBEStep( V_t )
    t = t - dt;
    
    if volatility_model > 1
        CalculateABCVectorsBE();
    end
    
    V_out_Smin = 0;
    if ~is_call
        if is_American_style
            V_out_Smin = X;
        else
            V_out_Smin = X * exp( -r*(T-t));
        end
    end
    
    RHS = V_t;
    RHS(i_first) = RHS(i_first) + A(i_first) * V_out_Smin;
    
    SolveTridiagonal( V_t );
    
    V_out(1) = V_out_Smin;
    
    if bc_Smax_Dirichlet
        V_out(n_price_levels) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = ApplyDividendAbsolute( V_in, D )
    if D < 0.001 * dS
        V_out = V_in;
        return;
    end
    
    iD_double = D / dS;
    iD_int = floor(iD_double);
    w0 = iD_double - iD_int;
    w1 = 1.0 - w0;
    
    if iD_int > n_price_levels
        iD_int = n_price_levels;
    end
    
    V_out(1 : iD_int + 1) = V_in(1);
    
    i = iD_int + 2;
    iold = 1;
    while i <= n_price_levels 
        V_out(i) = w0 * V_in(iold) + w1 * V_in(iold + 1);
        i = i + 1;
        iold = iold + 1;
    end
    
    if is_American_style
        for i = n_price_levels : -1 : 1
            if V_out(i) < payoff(i)
                V_out(i) = payoff(i);
            end
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = ApplyDividendMinPercentageAbsolute( V_in, iD )
    amount = 0.0;
    percentage = 0.0;
    
    if iD <= nD_absolute
        amount = div_amounts(iD);
    end
    
    if iD <= nD_percentage
        percentage = div_percentages(iD);
    end
    
    if percentage < 1e-10
        ApplyDividendAbsolute( V_in, amount );
        return
    end
    
    for i = 1 : n_price_levels
        D = percentage * S(i);
        if amount > 0.0 && D > amount
            D = amount;
        end
        
        Sdropped = S(i) - D;
        if Sdropped < 0.0
            Sdropped = 0.0;
        end
        
        if Sdropped <= S(1)
            V_out(i) = V_in(1);
        else
            inew_double = (Sdropped - S(1)) / dS + 1;
            inew = floor( inew_double );
            if inew >= n_price_levels
                inew = n_price_levels - 1;
            end
            w1 = inew_double - inew;
            w0 = 1.0 - w1;
            V_out(i) = w0 * V_in(inew) + w1 * V_in(inew+1);
            
        end
            
    end
    
    if is_American_style
        for i = n_price_levels : -1 : 1
            if V_out(i) < payoff(i)
                V_out(i) = payoff(i);
            end
        end
    end
    
end

end % FDMUniformGrid
