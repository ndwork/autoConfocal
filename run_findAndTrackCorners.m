
function run_findAndTrackCorners
  clear; close all; rng(1);

  img1 = rgb2gray( double(imread( './data/img1.jpg' ) ) / 255. );
  img2 = rgb2gray( double(imread( './data/img2.jpg' ) ) / 255. );

  N = 100;    % number of features
  buffer = 100;    % spacing between features
  searchWidth = 1000;
  [pts1,pts2] = findAndTrackCorners( img1, img2, 'N', N, ...
    'buffer', buffer, 'searchWidth', searchWidth );

  figure; imshow(img1,[]); labelImgPts( pts1 );  title('Image 1');
  figure; imshow(img2,[]); labelImgPts( pts2 );  title('Image 2');

  distThresh = 5;
  H12 = ransacDltHomographyFromPts2D( fliplr(pts1), fliplr(pts2), distThresh );

  aligned = projectImage( img1, H12 );
  figure; imshow( aligned, [] );
end

