% This software is offered under the GNU General Public License 3.0.  It 
% is offered without any warranty expressed or implied, including the 
% implied warranties of merchantability or fitness for a particular 
% purpose.

% Written by Nicholas Dwork (ndwork@stanford.edu) and 
% Gennifer Smith (gsmith9@stanford.edu).


function h = makeConfocalFunction( z, z0, zR )
  % z0 is the focal plane location
  % zR is the apparent Rayleigh Range

  tmp = (z-z0)/zR;
  h = 1 ./ ( tmp.^2 + 1 );

end
