
function [yShear,t,inlierIndxs] = ransacYShearAndTrans( pts1, pts2, thresh )
  % [yShear,t,inlierIndxs] = ransacYShearAndTrans( pts1, pts2, thresh )
  %
  % Uses RANSAC to find the yShear angle and vShift that aligns
  % pts2(:,1) = pts1(:,1) + t(1)
  % pts2(:,2) = tan(shear)*pts1(:,1) + pts2(:,2) + t(2);
  %
  % Inputs:
  % pts1,pts2 - 2D array with N rows and 2 columns
  %   M is the nubmer of points.  First/Second column is x/y location.
  % thresh - the RANSAC threshold
  %
  % Outputs:
  % yShear - shear angle in radians (to be used with shearImg)
  % t - the translation vector
  %
  % Written by Nicholas Dwork - Copyright 2016
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  [M, N] = size( pts1 );
  p = 0.99;

  sample_count = 0;
  bestNInliers = 0;
  nRansac = 1500;
  while nRansac > sample_count
    %disp(['Working on ', num2str(sample_count), ' of ', num2str(nRansac)]);
    indxs = randsample(M,2);  % 2 points determine yShear and t
    subset1 = pts1(indxs,:);
    subset2 = pts2(indxs,:);
    [yShear,t] = findYShearAndTrans( subset1, subset2 );

    aligned1 = pts1;
    aligned1(:,1) = pts1(:,1) + t(1);
    aligned1(:,2) = tan(yShear) * pts1(:,1) + pts1(:,2) + t(2);

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

  [yShear,t] = findYShearAndTrans( pts1(inlierIndxs,:), pts2(inlierIndxs,:) );
end


function [yShear,t] = findYShearAndTrans( pts1, pts2 )
  b = [ pts2(:,2) - pts1(:,2); pts2(:,1) - pts1(:,1); ];
  M = size(pts1,1);
  A = zeros(2*M,3);
  A(1:M,1) = pts1(:,1);
  A(1:M,3) = 1;
  A(M+1:end,2) = 1;
  v = A \ b;

  yShear = atan( v(1) );
  t = v(2:3);
end

