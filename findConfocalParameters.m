
function [z0, zR] = findConfocalParameters( bscan1, bscan2, dx, dz, shear, ...
  varargin )
  % [z0, zR] = findConfocalParameters( img1_dB, img2_dB, dx, dz, shear [, ...
  %    mask ] )
  %
  % Determines the confocal function parameters assuming that the only
  % difference between img1_dB and img2_dB is that the confocal function
  % has been vertically shifted.
  %
  % Inputs:
  % img1_dB - 2D array specifying img1 in decibels
  % img2_dB - 2D array specifying img2 in decibels aligned with img1.  It
  %   is assumed that all data in img1 overlaps with all data in img2
  % dx - (default 0) horizontal translation
  % dz - (default 0) vertical translation
  % shear - (default 0) amount to shear image 2
  %
  % Outputs:
  % z0 - the focal plane depth (in pixels)
  % zR - the Rayleigh range (in pixels)
  %
  % Written by Nicholas Dwork (ndwork@stanford.edu) - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  goodThresh = 8;
  noiseFraction = 0.05;

  defaultMask = ones( size(bscan1) );
  p = inputParser;
  p.addOptional( 'mask', defaultMask );
  p.parse( varargin{:} );
  mask = p.Results.mask;

  noiseIndx = floor( noiseFraction*size(bscan1,1) );
  tmp = bscan1( 1:noiseIndx, : );
  noiseMean1 = median( tmp(:) );
  deBiased1 = max( bscan1 - noiseMean1, min(bscan1(:)) );
  tmp = bscan2( 1:noiseIndx, : );
  noiseMean2 = median( tmp(:) );
  deBiased2 = max( bscan2 - noiseMean2, min(bscan2(:)) );
  
  bscan1_dB = intensity2dB( deBiased1 );
  bscan2_dB = intensity2dB( deBiased2 );

  bscan2_dB = project2onto1( bscan2_dB, dx, dz, shear, 'linear' );
  mask = project2onto1( mask, dx, dz, shear, 'nearest' );

  xs = 1:size(bscan1,2);
  vShifts = zeros(1,size(bscan1,2));
  rdx = round(dx);
  if rdx>=0
    vShifts(1:end-rdx) = xs(rdx+1:end)*tan(shear) + round(dz);
  else
    vShifts(-rdx+1:end) = xs(1:end+rdx)*tan(shear) + round(dz);
  end
  %vShifts = -vShifts;
  %vShifts is the amount to focal plane is shifted up in every column of 2

  gWidth = 19;
  gFilter = fspecial( 'gaussian', gWidth, 5 );
  smooth1 = imfilter( bscan1_dB, gFilter );
  smooth2 = imfilter( bscan2_dB, gFilter );
  se = strel( 'square', gWidth );
  mask = imerode( mask, se );

  smooth1 = smooth1(:,gWidth:end-gWidth);
  smooth2 = smooth2(:,gWidth:end-gWidth);
  mask = mask(:,gWidth:end-gWidth);
  mask(1:gWidth,:) = 0;
  vShifts = vShifts(gWidth:end-gWidth);

% LCol = 250;
% RCol = 270;
% smooth1 = smooth1(:,LCol:RCol);
% smooth2 = smooth2(:,LCol:RCol);
% mask = mask(:,LCol:RCol);
% vShifts = vShifts(LCol:RCol);

  diff_dB = smooth1 - smooth2;

  [nY, nX] = size( smooth1 );
  zImg1 = (1:nY)' * ones(1,nX);
  zImg2 = zImg1 + ones(nY,1) * vShifts;
  
  %normWeights = abs( smooth1 - smooth2 );
  minImg = min( smooth1, smooth2 );
  meanImg = (smooth1+smooth2)/2;
  %normWeights = minImg;
  normWeights = meanImg;
  %normWeights = diff_dB;
  normWeights( minImg < goodThresh | mask==0 ) = 0;

figure;  subplot(2,1,1);
c = mean(diff_dB,2);  plot( c(1:357), 'k', 'LineWidth', 2 )
subplot(2,1,2);  plot( normWeights(1:357), 'b', 'LineWidth', 2 );

  z0s = 2:nY-1;
  zRs = 1:floor(0.5*nY);
  nZ0 = numel(z0s);
  nZR = numel(zRs);
  costs = 99d999 * ones( nZ0, nZR );
  parfor j = 1:nZ0
    thisZ0 = z0s(j);
    if mod(j,50)==0,
      disp(['Working on ', num2str(j), ' of ', num2str(nZ0)]);
    end

    for i = 1:nZR
      thisZR = zRs(i);

      %model1 = intensity2dB( ( (zImg1-thisZ0) / thisZR ).^2 + 1 );
      %model2 = intensity2dB( ( (zImg2-thisZ0) / thisZR ).^2 + 1 );
      model1 = 10*log10( ( (zImg1-thisZ0) / thisZR ).^2 + 1 );
      model2 = 10*log10( ( (zImg2-thisZ0) / thisZR ).^2 + 1 );
      model = model2 - model1;
      costs(j,i) = norm( normWeights(:) .* (diff_dB(:) - model(:)), 1 );
    end
  end

  [~,minCostIndx] = min( costs(:) );
  [minZ0Indx, minZRIndx] = ind2sub( size(costs), minCostIndx );
  z0 = z0s( minZ0Indx );
  zR = zRs( minZRIndx );

showFeaturesOnImg( [minZRIndx, minZ0Indx], costs, [], 'scale', 2, 'color', 'y' );
model1 = 10*log10( ( (zImg1-z0) / zR ).^2 + 1 );
model2 = 10*log10( ( (zImg2-z0) / zR ).^2 + 1 );
model = model2 - model1;

figure; hold on;
c = mean(diff_dB,2);
plot(c,'k','LineWidth',2);
midCol = round( size( diff_dB, 2 ) / 2 );
dz_mm = 2.57 / 512;
trueZ0 = 1.2/dz_mm;
trueZR = 0.27/dz_mm;
trueModel1 = 10*log10( ( (zImg1-trueZ0) / trueZR ).^2 + 1 );
trueModel2 = 10*log10( ( (zImg2-trueZ0) / trueZR ).^2 + 1 );
trueModel = trueModel2 - trueModel1;
plot( trueModel(:,midCol), 'b', 'LineWidth', 2 );
plot( model(:,midCol), 'r', 'LineWidth', 2 );
disp(['z0(mm): ', num2str(z0*dz_mm)]);
disp(['zR(mm): ', num2str(zR*dz_mm)]);
plot( model(:,midCol), 'r', 'LineWidth', 2 );
legend('data','answer','model');
figure; plot( normWeights(:,midCol), 'LineWidth', 2 );
end


function proj2 = project2onto1( img2, dx, dz, shear, interpType )
  if nargin < 5
    interpType = 'linear';
  end

  sBScan = size(img2);
  [xs,ys] = meshgrid( 1:sBScan(2), 1:sBScan(1) );

  projCoords = zeros( 2, numel(xs) );
  projCoords(1,:) = xs(:) + dx;
  projCoords(2,:) = tan(shear)*xs(:) + ys(:) + dz;

  proj2 = interp2( xs, ys, img2, projCoords(1,:), projCoords(2,:), ...
    interpType, 0 );
  proj2 = reshape( proj2, sBScan );
end


