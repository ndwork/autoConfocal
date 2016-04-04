% This software is offered under the GNU General Public License 3.0.  It 
% is offered without any warranty expressed or implied, including the 
% implied warranties of merchantability or fitness for a particular 
% purpose.
%
% Written by Nicholas Dwork (ndwork@stanford.edu) and 
% Gennifer Smith (gsmith9@stanford.edu).


function [z0,zR,scaling] = findFocalPlane( img1, img2 )
  % img1 is a 2D array where intensity values are in decibels
  % img2 is a 2D array where intensity values are in decibels
  % z0 and zR are returned in pixels
  % scaling has no units
  
  goodThresh = median(img1(:)) + std(img1(:));

  vShift = findVertShift( img1, img2 );

  gWidth = 19;
  gFilter = fspecial( 'gaussian', gWidth, 5 );
  smooth1 = imfilter( img1, gFilter );
  smooth2 = imfilter( img2, gFilter );

  smooth1 = smooth1( gWidth:end-gWidth, gWidth:end-gWidth );
  smooth2 = smooth2( gWidth:end-gWidth, gWidth:end-gWidth );

  % Cut off edges due to corruptions
  smooth1 = smooth1( :, 50-gWidth: end-50+gWidth );
  smooth2 = smooth2( :, 50-gWidth: end-50+gWidth );

  sub1 = smooth1(1:end-vShift+1,:);
  sub2 = smooth2(vShift:end,:); % shifted up

  goodIndxs = find( sub1 > goodThresh );
  
  subDiff = sub1 - sub2;
  subDiffVec = subDiff(:);

  [nY,nX] = size( sub1 );
  zImg = (1:nY)' * ones(1,nX);

  z0s = 2:nY-1;
  zRs = 1:floor(0.5*nY);
  nZ0 = numel(z0s);
  nZR = numel(zRs);
  costs = 99d999 * ones( nZ0, nZR );
  parfor j = 1:nZ0
    thisZ0 = z0s(j);
    if mod(j,5)==0,
      disp(['Working on ', num2str(j), ' of ', num2str(nZ0)]);
    end

    for i = 1:nZR
      thisZR = zRs(i);

      h1 = 1 ./ ( ( (zImg-thisZ0) / thisZR ).^2 + 1 );
      h2 = 1 ./ ( ( (zImg+vShift-thisZ0) / thisZR ).^2 + 1 );
      model1 = 10*log10( h1 );
      model2 = 10*log10( h2 );
      model = model1 - model2;
      %model1 = 10*log10( ( (zImg+vShift-thisZ0) / thisZR ).^2 + 1 );
      %model2 = 10*log10( ( (zImg-thisZ0) / thisZR ).^2 + 1 );
      %model = thisScaling * ( model1 - model2 );
      %costs(j,i) = norm( subDiffVec(goodIndxs) - model(:), 1 );
      costs(j,i) = norm( subDiffVec(goodIndxs) - model(goodIndxs), 1 );

    end
  end

  [~,minCostIndx] = min( costs(:) );
  [minZ0Indx, minZRIndx] = ind2sub( size(costs), minCostIndx );
  z0 = z0s( minZ0Indx );
  zR = zRs( minZRIndx );
end

%thisZ0 = 512 / 2.57 * 1.4;  thisZR = 512 / 2.57 * 2*1.4*0.1;
%thisZR = 512 / 2.57 * 2*1.4*0.1;



function vShift = findVertShift( img1, img2 )
  nY = size( img1, 1 );
  nTest = floor(nY/2);

  costs = zeros(nTest,1);
  for i=1:nTest
    sub1 = img1(1:end-i+1,:);
    sub2 = img2(i:end,:);
    costs(i) = norm( sub1(:) - sub2(:), 2 ) / norm( sub1(:), 2 );
    
    
    %shifted = shiftImg( img2, [-i 0] );
    %costs(i) = norm( img1(:) - shifted(:), 2 );
  end

  [~,vShift] = min( costs );
  vShift = vShift - 1;    % Account for one based indexing in Matlab
end

