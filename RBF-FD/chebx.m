function x=chebx(rg,N);
  % Normal right cheb interval has length L=1 and starts at a=-1
  L=1; a=-1;  
  theta=(pi:-pi/2/(N-1):pi/2)';
  x = (rg(2)-rg(1))/L*(rg(1)-a+cos(theta)); 
end
