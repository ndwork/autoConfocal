
function [z0,zR,overlap] = findConfocal_yShearAndTrans( bscan1, bscan2, ...
  trueZ0_mm, trueZR_mm, dx_mm, dz_mm )
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
  gSig = 5;
  gWidth = 19;
  searchPercentage = 0.9;
  ransacThresh = 5;  % pixels

  bscan1_dB = intensity2dB( bscan1 );
  bscan2_dB = intensity2dB( bscan2 );

  h = fspecial('gaussian', gWidth, gSig);
  smooth1 = imfilter( bscan1_dB, h, 'symmetric', 'same' );
  smooth2 = imfilter( bscan2_dB, h, 'symmetric', 'same' );

  regWithPoints = 1;  % Choose registration algorithm

  if( regWithPoints )
    [pts1,pts2] = findAndTrackCorners( smooth1, smooth2, ...
     'w', 31, 'searchWidth', searchPercentage*size(smooth1) );
    figH1 = figure; imshowscale( bscan1_dB );  title( 'Features Image 1' );
    figH2 = figure; imshowscale( bscan2_dB );  title( 'Features Image 2' );

    sBScan = size( bscan1_dB );
    midX = ceil( (sBScan(2)+1)/2 );
    midY = ceil( (sBScan(1)+1)/2 );
    pts1(:,1) = pts1(:,1) - midX;  pts1(:,2) = pts1(:,2) - midY;
    pts2(:,1) = pts2(:,1) - midX;  pts2(:,2) = pts2(:,2) - midY;

    [yShear,t,inlierIndxs] = ransacYShearAndTrans( pts1, pts2, ransacThresh );

    pts2Plot1 = pts1;  pts2Plot2 = pts2;
    pts2Plot1(:,1) = pts2Plot1(:,1) + midX;  pts2Plot2(:,1) = pts2Plot2(:,1) + midX;
    pts2Plot1(:,2) = pts2Plot1(:,2) + midY;  pts2Plot2(:,2) = pts2Plot2(:,2) + midY;
    figure( figH1 ); labelImgPts( pts2Plot1, 'inlierIndxs', inlierIndxs );
    figure( figH2 ); labelImgPts( pts2Plot2, 'inlierIndxs', inlierIndxs );

    hShift = t(1);  vShift = t(2);
  else

    [vShift, hShift, rotation] = findTransRotWithPCC( bscan1_dB, bscan2_dB, ...
      'dTheta', 1*pi/180 );
    yShear = atan( rotation );
  end

  [z0, zR, overlap] = findConfocalParameters( bscan1, bscan2, ...
    hShift, vShift, yShear, trueZ0_mm, trueZR_mm, dx_mm, dz_mm );
end





