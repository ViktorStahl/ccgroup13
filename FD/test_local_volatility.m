clear;

profiling = 0;

if profiling
    profile on;
end

parameter_set = 1;

is_call = 1; % 1 for call, 0 for put
is_American_style = 0; % 0 for European, 1 for American style
barrier_type = 0; % 0 for vanilla option, 1 for up-and-out barrier
T = 1.0; %time to maturity
K = 100.0; %strike price
if barrier_type == 0
    B = 0; %barrier level
else
    B = 1.25 * K;
end
r = 0.03; %interest rate
q = 0.0; %dividend yield
volatility_model = 2; % 1 means constant volatility
volatility = 0.15; % volatility value

if parameter_set == 2
    T = 0.25;
    r = 0.1;
    volatility = 0.01;
end

%dividends
div_datetimes   = []; %[0.5 * T];
div_amounts     = []; %[5];  %absolute dividend anounts
div_percentages = []; %[ 0.05  ]; %proportional = percentage dividends

%settings for price grid
n_price_steps_per_strike = 400; 
Smax_multiplier = 2.0;
stick_strike_to_node = 1;

%settings for time grid
n_time_steps_per_year = 50;
min_time_steps_per_interval = 20;
n_Rannacher_sub_steps = 4; 

[S U payoff] = FDMUniformGrid( is_call, is_American_style, ...
    T, K, barrier_type, B, ...
    r, q, volatility_model, volatility, ...
    div_datetimes, div_amounts, div_percentages, ... 
    n_price_steps_per_strike, Smax_multiplier, stick_strike_to_node, ...
    n_time_steps_per_year, min_time_steps_per_interval, n_Rannacher_sub_steps);

m = 4; % grid refinement
l = 1.5; % extra Smax multiplier
[Sfine Ufine payoff_fine] = FDMUniformGrid( is_call, is_American_style, ...
    T, K, barrier_type, B, ...
    r, q, volatility_model, volatility, ...
    div_datetimes, div_amounts, div_percentages, ... 
    m * n_price_steps_per_strike, l * Smax_multiplier, stick_strike_to_node, ...
    m * n_time_steps_per_year, m * min_time_steps_per_interval, 2 * n_Rannacher_sub_steps);

n_price_levels = length(S);
n_fine = length(Sfine);

i_Smax = round( (S(n_price_levels) - Sfine(1)) / (Sfine(2) - Sfine(1)) + 1 );
if abs( Sfine(i_Smax) - S(n_price_levels) ) > 0.001 * (Sfine(2) - Sfine(1))
    error('Non-matched nodes of basic and refined grids');
end
Sref = Sfine(1: m: i_Smax);
Uref = Ufine(1: m: i_Smax);

dS = S(2) - S(1);


error_U = U - Uref;

error_U_relative(n_price_levels) = 0;
for i = 1 : n_price_levels
    if U(i) > 1e-2 * K
        error_U_relative(i) = error_U(i) / U(i);
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
U_requested = U( i_found );

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
U_str = ['U = [' sprintf('%0.15g ', U_requested) ']'];

disp(info); 
disp(info2);
if length(div_datetimes) > 0 
    disp(info3);
    disp(info4);
end

disp(S_str);
disp(U_str);

%figures
S_left = 0.6 * K;
if barrier_type == 0
    S_right = 1.6 * K;
else
    S_right = 1.25 * K;
end
out_step = 1;

i_left = round( S_left - S(1) ) / dS + 1;
i_right = round( S_right - S(1) ) / dS + 1;
jump = round( out_step / dS );

S_out = S(i_left : jump : i_right);
payoff_out = payoff(i_left : jump : i_right);

U_out = U(i_left : jump : i_right);

Uref_out = Uref(i_left : jump : i_right);

error_U_out = error_U(i_left : jump : i_right);
error_U_relative_out = error_U_relative(i_left : jump : i_right);

f1 = figure(1);
%set(0,'CurrentFigure',f1);
set(f1, 'Name', 'U');
plot(S_out, U_out, S_out, Uref_out, S_out, payoff_out);

f5 = figure(5);
set(f5,'Name', 'error U');
plot(S_out, error_U_out);

f6 = figure(6);
set(f6,'Name', 'error U relative');
plot(S_out, error_U_relative_out);

if profiling
    profile off;
    profile viewer;
    p = profile('info');
end
