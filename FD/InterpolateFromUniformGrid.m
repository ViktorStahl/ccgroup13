function [V] = InterpolateFromUniformGrid(S, S_FDM, V_FDM)
j_max = length(S);
i_max = length(S_FDM);
V(j_max) = 0;

dS = S_FDM(2) - S_FDM(1);

for j = 1:j_max
    i_S = floor(S(j) / dS + 0.1) + 1;
    if i_S > i_max-1
        i_S = i_max-1;
    end
    w = (S(j) - S_FDM(i_S)) / dS;

    if abs(w) < 1e-6
        V(j) = V_FDM(i_S);
    else
        V(j) = (1.0 - w) * V_FDM(i_S) + w * V_FDM(i_S + 1);
    end
end

end