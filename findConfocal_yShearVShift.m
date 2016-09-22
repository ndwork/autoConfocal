% This software is offered under the GNU General Public License 3.0.  It 
% is offered without any warranty expressed or implied, including the 
% implied warranties of merchantability or fitness for a particular 
% purpose.

function [z0,zR] = findConfocal_yShearVShift( bscan1_dB, bscan2_dB )
  % [z0,zR] = findConfcoal( bscan1_dB, bscan2_dB )
  %
  % Inputs:
  % bscan1_dB - first b-scan in decibels
  % bscan2_dB - second b-scan in decibels
  %
  % Outputs:
  % z0 - the focal plane location (in pixels)
  % zR - the Rayleigh range (in pixels)
  %
  % Written by Nicholas Dwork - Copyright 2016

  bscan1 = 10.^(bscan1_dB/10);
  bscan2 = 10.^(bscan2_dB/10);
  surfaceLocs1 = findSurfaceLocs( bscan1 );
  surfaceLocs2 = findSurfaceLocs( bscan2 );

  nLocs = numel(surfaceLocs1);
  nX = size( bscan1, 2 );
  midX = ceil( (nX+1) / 2 );
  pts1 = zeros(nLocs,2);       pts2 = zeros(nLocs,2);
  pts1(:,1) = (1:nLocs) - midX;  pts1(:,2) = surfaceLocs1;
  pts2(:,1) = (1:nLocs) - midX;  pts2(:,2) = surfaceLocs2;
  [yShear,vShift] = ransacYShearVShift( pts1, pts2, 5 );
  vShift = round( vShift );

  shearVShifts = pts1(:,1)*tan(yShear);

  gWidth = 19;
  gFilter = fspecial( 'gaussian', gWidth, 5 );
  smooth1 = imfilter( bscan1_dB, gFilter );
  smooth2 = imfilter( bscan2_dB, gFilter );

  mask = ones( size( bscan1_dB ) );
  mask = shearImg( mask, -yShear );
  sheared2 = shearImg( smooth2, -yShear );

  sub1 = smooth1(1:end-vShift,:);
  sub2 = sheared2(vShift+1:end,:); % shifted up
  subMask = mask(vShift+1:end,:); % shifted up

  se = strel( 'square', 7 );
  subMask = imerode( subMask, se );

% figure; imshow( sub1, [] );  title('aligned1 dB');
% figure; imshow( sub2, [] );  title('aligned2 dB');
% figure; imshow( mask, [] );  title('mask');

  [z0, zR] = findConfocalParameters( sub1, sub2, subMask, shearVShifts );
end





