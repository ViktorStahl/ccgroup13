function [ U delta gamma vega ] = BSeuExact( is_call, S, K, r, q, sigma, T )
%Black-Scholes formula for prices and Greeks of European call and put options
%
% http://en.wikipedia.org/wiki/Black%E2%80%93Scholes_model
%
% function [ U delta gamma vega ] = BSeuExact( is_call, S, K, r, q, sigma, T )
%
% Input arguments:
% is_call: 1 for call option, 0 for put option
% S: spot price
% K: strike price
% r: annualized interest rate
% q: annualized dividend yield
% sigma: annualized constant volatility
% T: time to maturity in years
%
% Output arguments:
% U: option price
% delta: option delta
% gamma: option gamma
% vega: option vega


sqrtT = sqrt(T);
sigma_sqrtT = sigma * sqrtT;
d1 = ( log(S/K) + (r - q + 0.5 * sigma * sigma) * T ) / sigma_sqrtT;
d2 = d1 - sigma_sqrtT;

exp_minus_qT = exp(-q * T);
K_exp_minus_rT = K * exp(-r * T);

if is_call
    delta = exp_minus_qT * normcdf(d1);
    U = S .* delta - K_exp_minus_rT * normcdf(d2);
else
    delta = - exp_minus_qT * normcdf(-d1);
    U = K_exp_minus_rT * normcdf(-d2) + S .* delta;
end

PDF_d1 = normpdf(d1);
gamma = exp_minus_qT * PDF_d1 ./ S / sigma_sqrtT;

vega = S .* PDF_d1 * sqrtT * exp_minus_qT;

end

