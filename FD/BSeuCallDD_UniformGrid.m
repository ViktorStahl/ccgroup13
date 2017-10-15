function [U]=BSeuCallDD_UniformGrid(S,K,T,r,sig,D,alpha)
%BSeuCallU_FDM_UniformGrid Summary of this function goes here
%   Detailed explanation goes here

is_call = 1; % 1 for call, 0 for put
is_American_style = 0; % 0 for European, 1 for American style
barrier_type = 0; % 0 for vanilla option, 1 for up-and-out barrier

B = 0; %barrier level
q = 0.0; %dividend yield
volatility_model = 1; % 1 means constant volatility

%dividends
div_datetimes = [alpha*T]; 
div_amounts   = [];   %absolute dividend anounts
div_percentages = [D]; %proportional = percentage dividends

%settings for price grid
n_price_steps_per_strike = 300; 
Smax_multiplier = 2.0;
stick_strike_to_node = 1;

%settings for time grid
n_time_steps_per_year = 100;
min_time_steps_per_interval = 20;
n_Rannacher_sub_steps = 8; 

[S_FDM U_FDM payoff] = FDMUniformGrid( is_call, is_American_style, ...
    T, K, barrier_type, B, ...
    r, q, volatility_model, sig, ...
    div_datetimes, div_amounts, div_percentages, ... 
    n_price_steps_per_strike, Smax_multiplier, stick_strike_to_node, ...
    n_time_steps_per_year, min_time_steps_per_interval, n_Rannacher_sub_steps);

U = InterpolateFromUniformGrid(S, S_FDM, U_FDM);
    
end

