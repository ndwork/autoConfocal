% This software is offered under the GNU General Public License 3.0.  It 
% is offered without any warranty expressed or implied, including the 
% implied warranties of merchantability or fitness for a particular 
% purpose.

% Written by Nicholas Dwork (ndwork@stanford.edu) and 
% Gennifer Smith (gsmith9@stanford.edu).


function f = makeFalloffFunction( z, lambda, deltaLambda, dLambda )
  % z = depth in mm
  % lambda = central frequency of light source
  % deltaLambda = wavelength per pixel in meters
  % dLambda = spectral resolution of spectrometer in meters

  zRD_mm = lambda*lambda/(4*deltaLambda) * 1000;
    % multiply by 1000 to convert from meters to mm

  zeta = pi/2 * z / zRD_mm;
  zetaSq = zeta .* zeta;

  w = dLambda / deltaLambda;

  sincZeta = sinc( zeta/pi );
  sincZetaSq = sincZeta .* sincZeta;
  f = sincZetaSq .* exp( -w*w./(2*log(2)) .* zetaSq );
end
