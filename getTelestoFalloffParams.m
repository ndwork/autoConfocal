
function [lambda0,deltaLambda,dLambda] = getTelestoFalloffParams()
  % [lambda0,deltaLambda,dLambda] = getTelestoFalloffParams()
  %
  % Outputs:
  % lambda0 - central frequency of light source
  % deltaLambda - wavelength per pixel in meters
  % dLambda - spectral resolution of spectrometer in meters
  %
  % Written by Nicholas Dwork (www.nicholasdwork.com) - Copyright 2018
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  lambda0 = 1325d-9;
  deltaLambda = 7.1875d-11;
  dLambda = 1.1055d-10;
end

