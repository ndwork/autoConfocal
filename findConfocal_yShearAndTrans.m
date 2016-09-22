
function [z0,zR] = findConfocal_yShearAndTrans( bscan1, bscan2 )
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

  noiseIndx = floor( noiseFraction*size(bscan1,1) );
  tmp = bscan1( noiseIndx:end, : );
  noiseMean = median( tmp(:) );

  bscan1_dB = 20*log10( bscan1 - noiseMean );
  bscan2_dB = 20*log10( bscan2 - noiseMean );
  
  h = fspecial('gaussian', gWidth, gSmooth);
  smooth1 = imfilter( bscan1_dB, h );
  smooth2 = imfilter( bscan2_dB, h );

  [pts1,pts2] = findAndTrackCorners( smooth1, smooth2, ...
    'w', 31, 'searchWidth', searchPercentage*size(smooth1) );
  figH1 = figure; imshow( bscan1_dB, [] ); labelImgPts( pts1 );
  figure; imshow( bscan2_dB, [] ); labelImgPts( pts2 );

  [yShear,t,inlierIndxs] = ransacYShearAndTrans( pts1, pts2, 10 );
  figure( figH1 );
  labelImgPts( pts1, 'inlierIndxs', inlierIndxs );

  % project image 2 onto image 1
  sBScan = size(bscan1_dB);
  [xs,ys] = meshgrid( 1:sBScan(2), 1:sBScan(1) );

  projCoords = zeros( 2, numel(xs) );
  projCoords(1,:) = xs(:) + t(1);
  projCoords(2,:) = tan(yShear)*xs(:) + ys(:) + t(2);

  proj2 = interp2( xs, ys, smooth2, projCoords(1,:), projCoords(2,:), ...
    'linear', 0 );
  proj2 = reshape( proj2, sBScan );
figure; imshow( smooth1, [] );
figure; imshow( proj2, [] );

  mask = ones( size( bscan1_dB ) );
  projMask = interp2( xs, ys, mask, projCoords(1,:), projCoords(2,:), ...
    'nearest', 0 );
  projMask = reshape( projMask, sBScan );
  projMask = imerode( projMask, strel( 'square', 7 ) );

% figure; imshow( sub1, [] );  title('aligned1 dB');
% figure; imshow( sub2, [] );  title('aligned2 dB');
% figure; imshow( mask, [] );  title('mask');

  xs = 1:sBScan(2);
  vShifts = zeros(1,sBScan(2));
  rt = round(t);
  if t(1)>=0
    vShifts(1:end-rt(1)) = xs(rt(1)+1:end)*tan(yShear) + t(2);
  else
    vShifts(-rt(1)+1:end) = xs(1:end+rt(1))*tan(yShear) + t(2);
  end
  vShifts = -vShifts;

  [z0, zR] = findConfocalParameters( smooth1, proj2, vShifts, 'mask', projMask );
end





