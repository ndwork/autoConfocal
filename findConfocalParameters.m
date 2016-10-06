
function [z0, zR] = findConfocalParameters( bscan1, bscan2, dx, dz, shear, ...
  trueZ0_mm, trueZR_mm )
  % [z0, zR] = findConfocalParameters( img1_dB, img2_dB, dx, dz, shear )
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

  %goodThresh = -90;
  noiseFraction = 0.05;
  %gWidth = 19;
  gWidth = [19 9];

  noiseIndx = floor( noiseFraction*size(bscan1,1) );
  tmp = bscan1( 1:noiseIndx, : );
  noiseMean1 = median( tmp(:) );
  deBiased1 = max( bscan1 - noiseMean1, min(bscan1(:)) );
  tmp = bscan2( 1:noiseIndx, : );
  noiseMean2 = median( tmp(:) );
  deBiased2 = max( bscan2 - noiseMean2, min(bscan2(:)) );

  bscan1_dB = intensity2dB( deBiased1 );
  bscan2_dB = intensity2dB( deBiased2 );

  surfaceLocs1 = findSurfaceLocs( deBiased1 );
surfaceFeatures = [ 1:size(deBiased1,2); surfaceLocs1'; ]';
showFeaturesOnImg( surfaceFeatures, bscan1_dB, [] );

  gFilter = fspecial( 'gaussian', gWidth, 5 );
  smooth1 = imfilter( bscan1_dB, gFilter, 'replicate' );
  smooth2 = imfilter( bscan2_dB, gFilter, 'replicate' );

  smooth1 = smooth1(:,gWidth:end-gWidth);
  smooth2 = smooth2(:,gWidth:end-gWidth);
  surfaceLocs1 = surfaceLocs1(gWidth:end-gWidth);
  mask = ones( size(smooth1) );
  mask(1:gWidth,:) = 0;
  mask(end-gWidth:end,:) = 0;

  smooth2 = project2onto1( smooth2, dx, dz, shear, 'linear' );
  mask = project2onto1( mask, dx, dz, shear, 'nearest' );
  mask(1:gWidth,:) = 0;
  for i=1:numel(surfaceLocs1)
    mask(1:surfaceLocs1(i),i) = 0;
  end
  mask = imerode( mask, strel( 'square', 19 ) );

figure; imshow( smooth1, [] );
figure; imshow( smooth2, [min(bscan1_dB(:)) max(bscan1_dB(:))] );

  xs = 1:size(bscan1,2);
  vShifts = zeros(1,size(bscan1,2));
  %vShifts is the amount the focal plane is shifted up in every column of 2
  rdx = round(dx);
  if rdx>=0
    vShifts(1:end-rdx) = xs(rdx+1:end)*shear + round(dz);
  else
    vShifts(-rdx+1:end) = xs(1:end+rdx)*shear + round(dz);
  end
  vShifts = vShifts(gWidth:end-gWidth);

xCoord = 270;
%xCoord = 400;
%xCoord = 30;
xWidth = 10;
figure; plotnice( mean(smooth1(:,xCoord-xWidth/2:xCoord+xWidth/2),2) )
hold on; plot( mean(smooth2(:,xCoord-xWidth/2:xCoord+xWidth/2),2), 'r', 'LineWidth', 2 )
disp(['vShifts(xCoord): ', num2str(vShifts(xCoord))]);
legend('bscan1', 'bscan2');

  diff_dB = smooth1 - smooth2;
figure; imshow( diff_dB .* mask, [] )

  [nY, nX] = size( smooth1 );
  z1s = (1:nY)' * ones(1,nX);
  z2s = z1s + ones(nY,1) * vShifts;

  %normWeights = abs( smooth1 - smooth2 );
  minImg = min( smooth1, smooth2 );
  meanImg = (smooth1+smooth2)/2;
  %normWeights = minImg > 10;
  normWeights = max(minImg,0);
normWeights = normWeights .* normWeights;
  %normWeights = sqrt( abs(max(meanImg,0)) );
  %normWeights = abs(diff_dB);
  %normWeights( minImg < goodThresh | mask==0 ) = 0;
  normWeights( mask==0 ) = 0;

sliver = 0;
if sliver == 1
  diff_dB = diff_dB(:,xCoord-xWidth/2:xCoord+xWidth/2);
  mask = mask(:,xCoord-xWidth/2:xCoord+xWidth/2);
  normWeights = normWeights(:,xCoord-xWidth/2:xCoord+xWidth/2);
  z1s = z1s(:,xCoord-xWidth/2:xCoord+xWidth/2);
  z2s = z2s(:,xCoord-xWidth/2:xCoord+xWidth/2);
end

  function cost = objective(x)
    thisZ0 = x(1);
    thisZR = x(2);
    model1 = 10*log10( ( (z1s-thisZ0) / thisZR ).^2 + 1 );
    %model1 = intensity2dB( ( (zImg1-thisZ0) / thisZR ).^2 + 1 );
    model2 = 10*log10( ( (z2s-thisZ0) / thisZR ).^2 + 1 );
    %model2 = intensity2dB( ( (zImg2-thisZ0) / thisZR ).^2 + 1 );
    model = model2 - model1;
    cost = norm( normWeights(:) .* (diff_dB(:) - model(:)), 1 );
  end

  exhaustive=0;
  if exhaustive==1
    z0s = 2:nY-1;
    zRs = 1:floor(0.5*nY);
    nZ0 = numel(z0s);
    nZR = numel(zRs);
    costs = 99d999 * ones( nZ0, nZR );
    %parfor j = 1:nZ0
for j = 199
      thisZ0 = z0s(j);
      if mod(j,50)==0,
        disp(['Working on ', num2str(j), ' of ', num2str(nZ0)]);
      end

      %for i = 1:nZR
      for i=58
        thisZR = zRs(i);

        %model1 = intensity2dB( ( (zImg1-thisZ0) / thisZR ).^2 + 1 );
        %model2 = intensity2dB( ( (zImg2-thisZ0) / thisZR ).^2 + 1 );
        model1 = 10*log10( ( (z1s-thisZ0) / thisZR ).^2 + 1 );
        model2 = 10*log10( ( (z2s-thisZ0) / thisZR ).^2 + 1 );
        model = model2 - model1;
%model = model(:,xCoord-xWidth/2:xCoord+xWidth/2);

figure; plotnice( mean(model,2) )
hold on; plot( mean(diff_dB,2), 'r', 'LineWidth', 2 )

        costs(j,i) = norm( normWeights(:) .* (diff_dB(:) - model(:)), 1 );
      end
    end

    [~,minCostIndx] = min( costs(:) );
    [minZ0Indx, minZRIndx] = ind2sub( size(costs), minCostIndx );
    z0 = z0s( minZ0Indx );
    zR = zRs( minZRIndx );

    showFeaturesOnImg( [minZRIndx, minZ0Indx], costs, [], 'scale', 2, 'color', 'y' );
  else

    x0 = [ floor(0.5*nY) floor(0.25*nY) ];  % x = [z0 zR]
    lb = [ 2 1 ];
    ub = [ nY-1 floor(0.5*nY) ];
    x = fmincon(@objective,x0,[],[],[],[],lb,ub);
    z0 = x(1);
    zR = x(2);

  end

model1 = 10*log10( ( (z1s-z0) / zR ).^2 + 1 );
model2 = 10*log10( ( (z2s-z0) / zR ).^2 + 1 );
model = model2 - model1;

figure; hold on;
c = mean(diff_dB.*mask,2);
plot(c,'k','LineWidth',2);
midCol = round( size( diff_dB, 2 ) / 2 );
dz_mm = 2.57 / 512;
trueZ0 = trueZ0_mm/dz_mm;
trueZR = trueZR_mm/dz_mm;
trueModel1 = 10*log10( ( (z1s-trueZ0) / trueZR ).^2 + 1 );
trueModel2 = 10*log10( ( (z2s-trueZ0) / trueZR ).^2 + 1 );
trueModel = trueModel2 - trueModel1;
plot( trueModel(:,midCol), 'b', 'LineWidth', 2 );
plot( model(:,midCol), 'r', 'LineWidth', 2 );
disp(['z0(mm): ', num2str(z0*dz_mm)]);
disp(['apparent zR(mm): ', num2str(zR*dz_mm)]);
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
  projCoords(2,:) = shear*xs(:) + ys(:) + dz;

  proj2 = interp2( xs, ys, img2, projCoords(1,:), projCoords(2,:), ...
    interpType, 0 );
  proj2 = reshape( proj2, sBScan );
end


