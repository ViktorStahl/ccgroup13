function [Vega] = BSeuCallVegaI_RBFFD(S,K,T,r,sig)
dsig = 1e-4;
Vl=BSeuCallUI_RBFFD(S,K,T,r,sig-0.5*dsig);
Vu=BSeuCallUI_RBFFD(S,K,T,r,sig+0.5*dsig);
Vega=(Vu-Vl)/dsig;
end