function sig = locvolI(s,t)
  
  sig = 0.15 + 0.15*(0.5+2*t).*(0.01*s-1.2).^2./((0.01*s).^2+1.44);