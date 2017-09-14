% This software is offered under the GNU General Public License 3.0.  It 
% is offered without any warranty expressed or implied, including the 
% implied warranties of merchantability or fitness for a particular 
% purpose.

% Written by Nicholas Dwork (ndwork@stanford.edu) and 
% Gennifer Smith (gsmith9@stanford.edu).


function mask = findNonZeroMus(image)
  % Input: image matrix
  % Output: binary image of the same size as the input, where the elements
  % where mu is 0 are set to 0
  
  dil = 10;
  erd = 25;
  
  [numR, numC] = size(image);

  mask = ones(numR, numC);

  for i = 1:numC
    skinLoc = findSurfaceLoc(image(:,i));
    mask(1:skinLoc, i) = 0;
  end
  
  mask = imdilate(mask, strel('square', dil));
  mask = imerode(mask, strel('square', erd));
  

end