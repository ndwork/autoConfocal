
function [yShear,vShift,inlierIndxs] = ransacYShearVShift( pts1, pts2, thresh )
  % [yShear,vShift] = ransacYShearVShift( pts1, pts2, thresh )
  %
  % Uses RANSAC to find the yShear angle and vShift that aligns
  % pts2 with pts1 such that pt2 = shear(pt1) + vShift;
  %
  % Inputs:
  % pts1,pts2 - 2D array with N rows and 2 columns
  %   N is the nubmer of points.  First/Second column is x/y location.
  % thresh - the RANSAC threshold
  %
  % Outputs:
  % yShear - shear angle in radians (to be used with shearImg)
  % vShift - vertical shift in pixels
  %
  % Written by Nicholas Dwork - Copyright 2016

  [M, N] = size( pts1 );
  p = 0.99;

  sample_count = 0;
  bestNInliers = 0;
  nRansac = 1500;
  while nRansac > sample_count
    %disp(['Working on ', num2str(sample_count), ' of ', num2str(nRansac)]);
    indxs = randsample(M,2);  % 2 points determine yShear and vShift
    subset1 = pts1(indxs,:);
    subset2 = pts2(indxs,:);
    [yShear,vShift] = findYShearVShift( subset1, subset2 );

    aligned1 = pts1;
    aligned1(:,2) = pts1(:,2) + pts1(:,1)*tan(yShear) + vShift;

    diffs = pts2 - aligned1;
    dists = sqrt( diffs(:,1).*diffs(:,1) + diffs(:,2).*diffs(:,2) );

    nInliers = sum( dists < thresh );
    if nInliers > bestNInliers
      bestNInliers = nInliers;
      inlierIndxs = find( dists < thresh );

      epsilon = 1 - nInliers / M;
      nRansac = log( 1 - p ) / log( 1 - (1-epsilon)^N );
    end

    sample_count = sample_count + 1;
  end

  [yShear,vShift] = findYShearVShift( ...
    pts1(inlierIndxs,:), pts2(inlierIndxs,:) );
end


function [yShear,vShift] = findYShearVShift( subset1, subset2 )
  nPts = size( subset1, 1 );
  A = [ subset1(:,1) ones(nPts,1) ];
  b = subset2(:,2) - subset1(:,2);
  vars = A \ b;
  yShear = atan( vars(1) );
  vShift = vars(2);
end

