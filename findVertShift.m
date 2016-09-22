
function vShift = findVertShift( img1, img2, mask )

  if nargin < 3, mask = ones( size(img1) ); end

  nY = size( img1, 1 );
  nTest = floor(nY/2);

  costs = zeros(nTest,1);
  for i=0:nTest-1
    sub1 = img1(1:end-i,:);
    sub2 = img2(i+1:end,:);
    subMask = mask(i+1:end,:);
    costs(i+1) = norm( sub1(subMask==1) - sub2(subMask==1), 2 ) / ...
      norm( sub1(subMask==1), 2 );
  end

  [~,vShift] = min( costs );
  vShift = vShift - 1;    % Account for one based indexing in Matlab
end

