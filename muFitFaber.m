
function muOut = muFitFaber( data, z, z0, zR, varargin )
  % muOut = muFitFaber( data, dz, z0, zR [, stdData] )
  % Based on the method described in "Quantitative measurement of
  % attenuation coefficients of weakly scattering media using optical coherence
  % tomography" by Faber et al., Optics Express, 2004
  %
  % Inputs:
  % data is a 1D array representing an A-scan OCT intensity values (perhaps averaged
  %   over a volume)
  % da is the vertical size of a pixel
  %
  % Optional Inputs:
  % stdData - array of size of data specifying standard deviation of each
  %   value of the A-scan
  %
  % Written by Nicholas Dwork - Copyright 2017

  defaultStdData = ones( size(data) );
  p = inputParser;
  p.addOptional( 'stdData', defaultStdData );
  p.parse( varargin{:} );
  stdData = p.Results.stdData;

  fitOptions = optimoptions( 'lsqnonlin', ...
    'Algorithm', 'levenberg-marquardt', ...
    'display', 'off', ...
    'TolFun', 1d-14, ...
    'TolX', 1d-14, ...
    'MaxFunEvals', 10000, ...
    'MaxIter', 10000 ...
  );

  h = makeConfocalFunction( z, z0, zR );

  function out = faberCost( x )
    % function is A * exp( -2*mu * z )
    A = x(1);
    mu = x(2);
    model = A * h .* exp( -2*mu * z );
    out = ( model(:) - data(:) ) ./ stdData(:);
  end

  x = lsqnonlin( @faberCost, [1 0], [], [], fitOptions );
  muOut= x(2);

  % Shift fit versus data
  model = x(1) * h .* exp( -2*muOut * z );
  figure; plotnice( z, model ); hold on; plotnice( z, data );
end

