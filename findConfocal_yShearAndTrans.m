
function [z0,zR] = findConfocal_yShearAndTrans( bscan1, bscan2, ...
  trueZ0_mm, trueZR_mm )
  % [z0,zR] = findConfocal_yShearAndTrans( bscan1, bscan2 )
  %
  % Inputs:
  % bscan1 - first b-scan
  % bscan2 - second b-scan
  %
  % Outputs:
  % z0 - the focal plane location (in pixels)
  % zR - the Rayleigh range (in pixels)
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  % Function parameters
  gSmooth = 5;
  gWidth = 19;
  searchPercentage = 0.8;
  noiseFraction = 0.9;
  ransacThresh = 15;  % pixels

  noiseIndx = floor( noiseFraction*size(bscan1,1) );
  tmp = bscan1( noiseIndx:end, : );
  noiseMean = median( tmp(:) );

  bscan1_dB = intensity2dB( bscan1 );
  bscan2_dB = intensity2dB( bscan2 );
  
  h = fspecial('gaussian', gWidth, gSmooth);
  smooth1 = imfilter( bscan1_dB, h, 'symmetric', 'same' );
  smooth2 = imfilter( bscan2_dB, h, 'symmetric', 'same' );

  [pts1,pts2] = findAndTrackCorners( smooth1, smooth2, ...
    'w', 31, 'searchWidth', searchPercentage*size(smooth1) );
  figH1 = figure; imshow( bscan1_dB, [] ); labelImgPts( pts1 );
  figure; imshow( bscan2_dB, [] ); labelImgPts( pts2 );

  [yShear,t,inlierIndxs] = ransacYShearAndTrans( pts1, pts2, ransacThresh );
  figure( figH1 );
  labelImgPts( pts1, 'inlierIndxs', inlierIndxs );

  [z0, zR] = findConfocalParameters( bscan1, bscan2, t(1), t(2), yShear, ...
    trueZ0_mm, trueZR_mm );
end





