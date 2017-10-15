function gamma = CalculateGamma(V,dS)
    n_price_levels = length(V);
    dS_squared = dS*dS;
    gamma(n_price_levels) = 0.0; 
    for i=2:n_price_levels-1
        gamma(i) = (V(i+1) - 2*V(i) + V(i-1)) / dS_squared;
    end
    gamma(1) = gamma(2);
    gamma(n_price_levels) = gamma(n_price_levels-1);

end
