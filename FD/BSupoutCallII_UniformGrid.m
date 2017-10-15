function [U]=BSupoutCallII_UniformGrid(S,K,T,r,sig,B)
%BSeuCallU_FDM_UniformGrid Summary of this function goes here
%   Detailed explanation goes here

is_call = 1; % 1 for call, 0 for put
is_American_style = 0; % 0 for European, 1 for American style
barrier_type = 1; % 0 for vanilla option, 1 for up-and-out barrier

q = 0.0; %dividend yield
volatility_model = 1; % 1 means constant volatility

%dividends
div_datetimes = []; 
div_amounts   = [];   %absolute dividend anounts
div_percentages = []; %proportional = percentage dividends

%settings for price grid
n_price_steps_per_strike = 40000; 
Smax_multiplier = 1.0; %no effect on barrier options
stick_strike_to_node = 1;

%settings for time grid
n_time_steps_per_year = 4800;
min_time_steps_per_interval = 1200;
n_Rannacher_sub_steps = 8; 

[S_FDM U_FDM payoff] = FDMUniformGrid( is_call, is_American_style, ...
    T, K, barrier_type, B, ...
    r, q, volatility_model, sig, ...
    div_datetimes, div_amounts, div_percentages, ... 
    n_price_steps_per_strike, Smax_multiplier, stick_strike_to_node, ...
    n_time_steps_per_year, min_time_steps_per_interval, n_Rannacher_sub_steps);

U = InterpolateFromUniformGrid(S, S_FDM, U_FDM);

j_max = length(S);
for j=1:j_max
    if S(j) >= B 
        U(j) = 0;
    end
end

   
end

