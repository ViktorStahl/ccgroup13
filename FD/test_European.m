clear; 

parameter_set = 2;

is_call = 1; % 1 for call, 0 for put
is_American_style = 0; % 0 for European, 1 for American style
barrier_type = 1; % 0 for vanilla option, 1 for up-and-out barrier
T = 1.0; %time to maturity
K = 100.0; %strike price
if barrier_type == 0
    B = 0; %barrier level
else
    B = 1.25 * K;
end
r = 0.03; %interest rate
q = 0.0; %dividend yield
volatility_model = 1; % 1 means constant volatility
volatility = 0.15; % volatility value

if parameter_set == 2
    T = 0.25;
    r = 0.1;
    volatility = 0.01;
end

%dividends
div_datetimes = []; %[0.5 * T];
div_amounts   = []; %[5]; %absolute dividend anounts
div_percentages = [];     %proportional = percentage dividends

%settings for price grid
n_price_steps_per_strike = 1000; 
Smax_multiplier = 1.2;
stick_strike_to_node = 1;

%settings for time grid
n_time_steps_per_year = 100;
min_time_steps_per_interval = 25;
n_Rannacher_sub_steps = 4; 

[S U_FDM payoff] = FDMUniformGrid( is_call, is_American_style, ...
    T, K, barrier_type, B, ...
    r, q, volatility_model, volatility, ...
    div_datetimes, div_amounts, div_percentages, ... 
    n_price_steps_per_strike, Smax_multiplier, stick_strike_to_node, ...
    n_time_steps_per_year, min_time_steps_per_interval, n_Rannacher_sub_steps);

dS = S(2) - S(1);
delta_FDM = CalculateDelta(U_FDM, dS);
gamma_FDM = CalculateGamma(U_FDM, dS);

volatility2 = volatility * 1.0001;
[S2 U_FDM2 payoff2] = FDMUniformGrid( is_call, is_American_style, ...
    T, K, barrier_type, B, ...
    r, q, volatility_model, volatility2, ...
    div_datetimes, div_amounts, div_percentages, ... 
    n_price_steps_per_strike, Smax_multiplier, stick_strike_to_node, ...
    n_time_steps_per_year, min_time_steps_per_interval, n_Rannacher_sub_steps);

vega_FDM = (U_FDM2 - U_FDM) / (volatility2 - volatility);

n_price_levels = length(S);

if barrier_type == 0
    [U delta gamma vega] = BSeuExact( is_call, S, K, r, q, volatility, T );
else
    U = OutBarrierExact(S, B, K, 0, T, volatility, r, q, is_call, 'U');
    delta = CalculateDelta(U, dS);
    gamma = CalculateGamma(U, dS);
    
    U2 = OutBarrierExact(S, B, K, 0, T, volatility2, r, q, is_call, 'U');
    vega = (U2 - U) / (volatility2 - volatility);
end


error_U = U_FDM - U;
error_delta = delta_FDM - delta;
error_gamma = gamma_FDM - gamma;
error_vega = vega_FDM - vega;

error_U_relative(n_price_levels) = 0;
error_delta_relative(n_price_levels) = 0;
error_gamma_relative(n_price_levels) = 0;
error_vega_relative(n_price_levels) = 0;
for i = 1 : n_price_levels
    if U(i) > 1e-2 * K
        error_U_relative(i) = error_U(i) / U(i);
        error_delta_relative(i) = error_delta(i) / delta(i);
        error_gamma_relative(i) = error_gamma(i) / gamma(i);
        error_vega_relative(i) = error_vega(i) / vega(i);
    end
end


%values
S_requested_vals = [90.0 100.0 110.0];
n_requested = length(S_requested_vals);
i_found(n_requested) = 0;
for j = 1 : n_requested
    i_found(j) = round( S_requested_vals(j) - S(1) ) / dS + 1;
end
S_requested = S( i_found );
U_FDM_requested = U_FDM( i_found );
delta_FDM_requested = delta_FDM( i_found );
gamma_FDM_requested = gamma_FDM( i_found );
vega_FDM_requested = vega_FDM( i_found );

error_U_requested = error_U_relative( i_found );
error_delta_requested = error_delta_relative( i_found );
error_gamma_requested = error_gamma_relative( i_found );
error_vega_requested = error_vega_relative( i_found );


if is_American_style
    info = 'American ';
else
    info = 'European ';
end

if barrier_type == 1
    info = [info 'up-and-out '];
end

if is_call
    info = [info 'call option'];
else
    info = [info 'put option'];
end

if length(div_datetimes) > 0 
    info = [info ' with dividends'];
end

info2 = sprintf('r=%g, q=%g, sigma=%g, K=%g, T=%g', r, q, volatility, K, T);
if length(div_datetimes) > 0 
    info3 = ['dividend datetimes: ' sprintf('%g ', div_datetimes)];
    info4 = ['dividend ratios: ' sprintf('%g ', div_percentages)];
end

S_str = ['S = [' sprintf('%0.15g ', S_requested) ']'];
U_str = ['U = [' sprintf('%0.15g ', U_FDM_requested) ']'];
delta_str = ['delta = [' sprintf('%0.15g ', delta_FDM_requested) ']'];
gamma_str = ['gamma = [' sprintf('%0.15g ', gamma_FDM_requested) ']'];
vega_str = ['vega = [' sprintf('%0.15g ', vega_FDM_requested) ']'];

error_U_str = ['error_U = [' sprintf('%0.15g ', error_U_requested) ']'];
error_delta_str = ['error_delta = [' sprintf('%0.15g ', error_delta_requested) ']'];
error_gamma_str = ['error_gamma = [' sprintf('%0.15g ', error_gamma_requested) ']'];
error_vega_str = ['error_vega = [' sprintf('%0.15g ', error_vega_requested) ']'];


disp(info); 
disp(info2);
if length(div_datetimes) > 0 
    disp(info3);
    disp(info4);
end

disp(S_str);
disp(U_str);
disp(delta_str);
disp(gamma_str);
disp(vega_str);

disp(error_U_str);
disp(error_delta_str);
disp(error_gamma_str);
disp(error_vega_str);


%figures
S_left = 0.6 * K;
if barrier_type == 0
    S_right = 1.6 * K;
else
    S_right = 1.25 * K;
end

if S_right > S(n_price_levels)
    S_right = S(n_price_levels);
end

out_step = 1;

i_left = round( S_left - S(1) ) / dS + 1;
i_right = round( S_right - S(1) ) / dS + 1;
jump = round( out_step / dS );

S_out = S(i_left : jump : i_right);
payoff_out = payoff(i_left : jump : i_right);

U_FDM_out = U_FDM(i_left : jump : i_right);
delta_FDM_out = delta_FDM(i_left : jump : i_right);
gamma_FDM_out = gamma_FDM(i_left : jump : i_right);
vega_FDM_out = vega_FDM(i_left : jump : i_right);

U_out = U(i_left : jump : i_right);
delta_out = delta(i_left : jump : i_right);
gamma_out = gamma(i_left : jump : i_right);
vega_out = vega(i_left : jump : i_right);

error_U_out = error_U(i_left : jump : i_right);
error_U_relative_out = error_U_relative(i_left : jump : i_right);

f1 = figure(1);
%set(0,'CurrentFigure',f1);
set(f1, 'Name', 'U');
plot(S_out, U_FDM_out, S_out, U_out, S_out, payoff_out);

f2 = figure(2);
set(f2, 'Name', 'delta');
plot(S_out, delta_FDM_out, S_out, delta_out);

f3 = figure(3);
set(f3, 'Name', 'gamma');
plot(S_out, gamma_FDM_out, S_out, gamma_out);
 
f4 = figure(4);
set(f4, 'Name', 'vega');
plot(S_out, vega_FDM_out, S_out, vega_out);
 
f5 = figure(5);
set(f5,'Name', 'error U');
plot(S_out, error_U_out);

f6 = figure(6);
set(f6,'Name', 'error U relative');
plot(S_out, error_U_relative_out);


