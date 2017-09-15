% This software is offered under the GNU General Public License 3.0.  It 
% is offered without any warranty expressed or implied, including the 
% implied warranties of merchantability or fitness for a particular 
% purpose.

% Written by Nicholas Dwork (ndwork@stanford.edu) and 
% Gennifer Smith (gsmith9@stanford.edu).

function skinLoc = findSurfaceLoc( line )

  n = numel(line);
  sortedLine =  sort(line);
  low = mean( sortedLine(1:round(0.1*n)) );
  high = median( sortedLine(round(0.9*n):end) );
  %high = mean( sortedLine(round(0.9*n):end) );

  thresh = ( high + low ) * 0.5;

  % Find the first point that is below the threshold and the first (after
  % the threshold) that is above the threshold  
  mask = line > thresh;
  maskDiffs = diff(mask);
  skinLoc = find( maskDiffs==1, 1 );
end
