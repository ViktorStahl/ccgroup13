function [phi,phi1,phi2]=iq(epsil,r,nprime,dim)
%
% NPRIME is a string defining which operator to use on the basis function
%
% DIM is the dimension for the partial derivative if nprime is '0','1'...,'4'
% and DIM is the number of space dimensions if nprime is 'L' or 'L2'
%
if nargin<=3
  dim=1;
end

if (nprime(1)=='L')
  if (size(r,3)==1)
    nd=dim;
  else
    nd=size(r,3)-1;
  end
end

phi=zeros(size(r,1),size(r,2));

tmp = 1+(epsil*sq(r(:,:,1))).^2;

if nprime(1)=='0'

  phi = 1./tmp;

elseif nprime(1)=='1'
 
  phi = -2*epsil^2*sq(r(:,:,dim+1))./tmp.^2;

elseif nprime=='2'

  phi = -2*epsil^2./tmp.^2 + 8*epsil^4*sq(r(:,:,dim+1)).^2./tmp.^3;

elseif nprime(1)=='3'

  phi = 24*epsil^4*sq(r(:,:,dim+1))./tmp.^3 - ...
        48*epsil^6*sq(r(:,:,dim+1)).^3./tmp.^4;

elseif nprime(1)=='4'

  phi = 24*epsil^4./tmp.^3 - ...
        288*epsil^6*sq(r(:,:,dim+1)).^2./tmp.^4 + ...
        384*epsil^8*sq(r(:,:,dim+1)).^4./tmp.^5;

elseif nprime(1)=='L' & length(nprime)==1

phi = -2*nd*epsil^2./tmp.^2 + 8*epsil^4*sq(r(:,:,1)).^2./tmp.^3;

elseif nprime(1:2)=='L2'

phi = 8*nd*(nd+2)*epsil^4./tmp.^3 - ...
      96*(nd+2)*epsil^6*sq(r(:,:,1)).^2./tmp.^4 + ...
                           384*epsil^8*sq(r(:,:,1)).^4./tmp.^5;
else
  error('Error in input argument nprime to function iq')
end 


function r=sq(r)
 r=squeeze(r);
