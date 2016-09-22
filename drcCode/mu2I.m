% This software is offered under the GNU General Public License 3.0.  It 
% is offered without any warranty expressed or implied, including the 
% implied warranties of merchantability or fitness for a particular 
% purpose.

% Written by Nicholas Dwork (ndwork@stanford.edu) and 
% Gennifer Smith (gsmith9@stanford.edu).


function I = mu2I( mu, z, z0, zR, alpha, beta, L0, ...
	lambda, deltaLambda, dLambda )
  % I = mu2I( mu, z, z0, zR, alpha, beta, L0, ...
	%   lambda, deltaLambda, dLambda );
  %
  % Inputs:
  % mu - 1D array representing vector of attenuation coefficients
  % z - 1D represents the depth at the front of each dz layer
  % z0 - scalar representing depth of focal length
  % zR - scalar representing apparent Rayleigh range
  %
  % Optional Inputs:
  % alpha - scalar representing back reflection (default is 1)
  % beta - scalar representing detector efficiency (default is 1)
  % L0 - scalar reprenseting initial light intensity (default is 1)
  % Note:  alpha, beta, L0 are all constants of proportionality
  % Note: lambda, deltaLambda, dLambda are all parameters for falloff
  %
  % Outputs:
  % I - represents the intensity in the middle of each dz layer

  dz = z(2:end) - z(1:end-1);
  dz = [ dz; dz(end) ];  % symmetric boundary condition

  [M,N] = size(mu);
  I = zeros( M, N );

  if nargin > 4
    k = alpha*beta*L0;
    if numel(k) == 0, k=1; end;
  else
    k = 1;
  end

  h = makeConfocalFunction( z, z0, zR );
  if nargin > 7
    f = makeFalloffFunction( z, lambda, deltaLambda, dLambda );
    hf = h .* f;
  else
    hf = h;
  end

  for j=1:N
    tmp = mu(1,j) * dz(1)/2;
    for i=1:M
      I(i,j) = k*mu(i,j) .* exp( -2 * tmp );
      if i < M
        tmp = tmp + mu(i,j)*dz(i)/2 + mu(i+1,j)*dz(i+1)/2;
      end
    end
    I(:,j) = dz .* I(:,j);  % Integration over a pixel
    I(:,j) = I(:,j) .* hf;
  end
end
