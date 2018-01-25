
function run_falloffSensitivity
  clear; close all; rng(1);
  addpath(genpath(pwd))

  %datacase =  1;  % Structured phantom with rotation and translation
  %datacase = 11;  % Zeiss retina
  datacase = 12;  % Rabbit eye
  %datacase = 14;  % Structured phantom 2 (of sensitivity analysis ) 
                   % with rotation and translation
  %datacase = 15;  % Sensitivity analysis phantom
  %datacase = 16;  % Just noise
  %datacase = 17;  % Focal plane high above sample
  [ bscan1, bscan2, dx_mm, dz_mm, noisePower, lambda0, deltaLambda, dLambda, ...
    trans,trueZ0_mm,trueZR_mm ] = loadDataCase( datacase );                                %#ok<ASGLU>

  % Remove existing falloff effect
  nRows = size( bscan1, 1 );
  nCols = size( bscan1, 2 );
  z = (0:nRows-1) / (nRows+1) * (2.57*nRows/512);
  f = makeFalloffFunction( z, lambda0, deltaLambda, dLambda );
  bscan1_noF = bscan1;  bscan2_noF = bscan2;
  for i=1:nCols
   bscan1_noF(:,i) = bscan1(:,i) ./ f';
   bscan2_noF(:,i) = bscan2(:,i) ./ f';
  end

  [~,~,~,vShift,hShift,yShear] = findConfocal_yShearAndTrans( bscan1, bscan2, ...
    lambda0, deltaLambda, dLambda, trueZ0_mm, trueZR_mm, dx_mm, dz_mm );
  close all;
%save( 'junk.mat', 'vShift', 'hShift', 'yShear' );
%load junk.mat;

  sigmas = 100:25:800;
  nSigmas = numel( sigmas );
  z0s = zeros( nSigmas, 1 );
  zRs = zeros( nSigmas, 1 );
  falloffRates = zeros( nSigmas, 1 );
  for sigmaIndx = 1:nSigmas
    disp([ 'Working on ', num2str(sigmaIndx), ' of ', num2str(nSigmas) ]);
    sigma = sigmas( sigmaIndx );

    zIndxs = 0:(nRows-1);
    newFalloff = exp( -zIndxs.^2 / (2*sigma*sigma) );
    newFalloff = newFalloff / max( newFalloff(:) );

    newBscan1 = zeros( size( bscan1 ) );
    newBscan2 = zeros( size( bscan2 ) );
    for i=1:nCols
      newBscan1(:,i) = bscan1_noF(:,i) .* newFalloff';
      newBscan2(:,i) = bscan2_noF(:,i) .* newFalloff';
    end

    [z0, zR, overlap_mmSq] = findConfocalParameters( newBscan1, newBscan2, ...
      hShift, vShift, yShear, lambda0, deltaLambda, dLambda, trueZ0_mm, trueZR_mm, ...
      dx_mm, dz_mm );                                                                      %#ok<ASGLU>

    z0s( sigmaIndx ) = z0;
    zRs( sigmaIndx ) = zR;

    % 6 dB is 50% of the amplitude
    z3dB = find( newFalloff < 0.708, 1, 'first' );
    if numel( z3dB ) > 0
      z3dB_mm = z3dB * dz_mm;
      falloffRates( sigmaIndx ) = 3 / z3dB_mm;  % decibels per mm
    else
      z1dB = find( newFalloff < 0.891, 1, 'first' );
      z1dB_mm = z1dB * dz_mm;
      falloffRates( sigmaIndx ) = 1 / z1dB_mm;  % decibels per mm
    end

    close all;
  end

  z0s_mm = z0s * dz_mm;
  zRs_mm = zRs * dz_mm;

  figure; plotnice( sigmas, z0s_mm );  title('z0s in mm');
  figure; plotnice( sigmas, zRs_mm );  title('zRs in mm');

  figure; plotnice( falloffRates, z0s_mm, 'k' );  title('z0s per falloff rate');
  figure; plotnice( falloffRates, zRs_mm, 'k' );  title('zRs per falloff rate');

end
