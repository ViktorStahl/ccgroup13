function [U]=BSeuLocVolII_UniformGrid(S,K,T,r,locvolfun)
%BSeuCallU_FDM_UniformGrid Summary of this function goes here
%   Detailed explanation goes here

is_call = 1; % 1 for call, 0 for put
is_American_style = 0; % 0 for European, 1 for American style
barrier_type = 0; % 0 for vanilla option, 1 for up-and-out barrier

B = 0; %barrier level
q = 0.0; %dividend yield
volatility_model = 3; % 1 means constant volatility

%dividends
div_datetimes = []; 
div_amounts   = [];   %absolute dividend anounts
div_percentages = []; %proportional = percentage dividends

%settings for price grid
n_price_steps_per_strike = 400; %1000; 
Smax_multiplier = 2.0;
stick_strike_to_node = 1;

%settings for time grid
n_time_steps_per_year = 100;
min_time_steps_per_interval = 25; 
n_Rannacher_sub_steps = 4; 

j_max = length(S);
U(j_max) = 0.0;
for j=1:j_max
    volatility_parameter = {locvolfun S(j)};
    [S_FDM U_FDM payoff] = FDMUniformGrid( is_call, is_American_style, ...
        T, K, barrier_type, B, ...
        r, q, volatility_model, volatility_parameter, ...
        div_datetimes, div_amounts, div_percentages, ... 
        n_price_steps_per_strike, Smax_multiplier, stick_strike_to_node, ...
        n_time_steps_per_year, min_time_steps_per_interval, n_Rannacher_sub_steps);

    dS = S_FDM(2) - S_FDM(1);
    i_S = floor(S(j) / dS + 0.1) + 1;
    i_max = length(S_FDM);
    if i_S > i_max-1
        i_S = i_max-1;
    end
    w = (S(j) - S_FDM(i_S))/dS;

    if abs(w) < 1e-6
        U(j) = U_FDM(i_S);
    else
        U(j) = (1.0 - w) * U_FDM(i_S) + w * U_FDM(i_S + 1);
    end
end

end

