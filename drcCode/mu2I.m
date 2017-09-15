% This software is offered under the GNU General Public License 3.0.  It 
% is offered without any warranty expressed or implied, including the 
% implied warranties of merchantability or fitness for a particular 
% purpose.

% Written by Nicholas Dwork (ndwork@stanford.edu) and 
% Gennifer Smith (gsmith9@stanford.edu).


function I = mu2I( mu, z, z0, zR, alpha, beta, L0 )
  % z represents the depth at the front of each dz layer
  % I represents the intensity in the middle of each dz layer

  dz = z(2:end) - z(1:end-1);
  dz = [ dz; dz(end) ];  % symmetric boundary condition

  I = zeros( numel(z), 1 );

  k = alpha*beta*L0;
  if numel(k) == 0, k=1; end;

  tmp = mu(1) * dz(1)/2;
  for i=1:numel(I)
    I(i) = k*mu(i) .* exp( -2 * tmp );
    if i < numel(I)
      tmp = tmp + mu(i)*dz(i)/2 + mu(i+1)*dz(i+1)/2;
    end
  end
  I = dz .* I;  % Integration over a pixel

  h = makeConfocalFunction( z, z0, zR );
  I = I .* h;
end
