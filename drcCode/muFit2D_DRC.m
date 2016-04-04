% This software is offered under the GNU General Public License 3.0.  It 
% is offered without any warranty expressed or implied, including the 
% implied warranties of merchantability or fitness for a particular 
% purpose.
%
% Written by Nicholas Dwork (ndwork@stanford.edu) and 
% Gennifer Smith (gsmith9@stanford.edu).


function muFit = muFit2D_DRC( I, z, z0, zR, noisePower, ...
  lambda, deltaLambda, dLambda )

  [M, N] = size( I );
  muFit = zeros( M, N );

  h = makeConfocalFunction( z, z0, zR );

  if nargin > 5
    f = makeFalloffFunction( z, lambda, deltaLambda, dLambda );
    hf = h .* f;
  else
    hf = h;
  end

  parfor j=1:N
    line = I(:,j);
    dz = z(2) - z(1);
    muFit(:,j) = muFitDRC( line, dz, hf, noisePower );
  end

end

