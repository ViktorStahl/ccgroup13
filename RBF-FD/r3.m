function [phi]=r3(epsil,r,nprime,dim);

if nargin<=3
  dim=1;
end

phi=zeros(size(r,1),size(r,2));

if nprime(1)=='0'
 
  phi= sq(r(:,:,1)).^3;

elseif nprime(1)=='1'

phi= 3*sq(r(:,:,1).*r(:,:,dim+1));

elseif nprime(1)=='2'

% Make sure that phi is computed correctly for r=0
  mask=(r(:,:,1)==0);  
  phi=3*(sq(r(:,:,1))+sq(r(:,:,dim+1).^2./(r(:,:,1)+mask)));

elseif nprime(1)=='L' & length(nprime)==1

nd = size(r,3)-1;
phi = 3*(nd+1)*sq(r(:,:,1));

else
  error('Error in input argument nprime to function r3')
end

function r=sq(r)
 r=squeeze(r);
