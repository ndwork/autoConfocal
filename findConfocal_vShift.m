
function [z0,zR] = findConfocal_vShift( bscan1, bscan2, trueZ0_mm, ...
  trueZR_mm, dx_mm, dz_mm )
  % [z0,zR] = findConfocal_vShift( bscan1, bscan2, trueZ0_mm, trueZR_mm, ...
  %   dx_mm, dz_mm )
  % Assumes img2 and img1 are related by a vertical translation
  %
  % Inputs:
  % bscan1 is a 2D array
  % bscan2 is a 2D array
  % mask is an optional input specifying which pixels in bscan2 are valid
  %
  % Outputs:
  % z0 and zR are returned in pixels
  %
  % Written by Nicholas Dwork (ndwork@stanford.edu)
  %
  % This software is offered under the GNU General Public License 3.0.  It 
  % is offered without any warranty expressed or implied, including the 
  % implied warranties of merchantability or fitness for a particular 
  % purpose.

  bscan1_dB = intensity2dB( bscan1 );
  bscan2_dB = intensity2dB( bscan2 );
  vShift = findVertShift( bscan1_dB, bscan2_dB );

  [z0, zR] = findConfocalParameters( bscan1, bscan2, 0, vShift, 0, ...
    trueZ0_mm, trueZR_mm, dx_mm, dz_mm );
end

%thisZ0 = 512 / 2.57 * 1.4;  thisZR = 512 / 2.57 * 2*1.4*0.1;
%thisZR = 512 / 2.57 * 2*1.4*0.1;


