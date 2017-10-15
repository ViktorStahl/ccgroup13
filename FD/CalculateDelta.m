function delta = CalculateDelta(V,dS)
    n_price_levels = length(V);
    dS_2 = dS * 2.0;
    delta(n_price_levels) = 0.0; 
    for i=2:n_price_levels-1
        delta(i) = (V(i+1) - V(i-1)) / dS_2;
    end
    delta(1) = delta(2);
    delta(n_price_levels) = delta(n_price_levels-1);
end
