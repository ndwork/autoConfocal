
function run_findConfocal
  clear; close all; rng(1);
  addpath(genpath(pwd))

  datacase = 8;
  [bscan1,bscan2,dz_mm,noisePower,trans,trueZ0_mm,trueZR_mm] = ...
    loadDataCase( datacase );

  switch trans
    case 'vShift'
      [z0,zR] = findConfocal_vShift( bscan1, bscan2 );
    case 'vShiftRot'
      [z0,zR] = findConfocal_vShiftRot( bscan1, bscan2 );
    case 'yShearAndTrans'
      [z0,zR] = findConfocal_yShearAndTrans( bscan1, bscan2, ...
        trueZ0_mm, trueZR_mm );
  end

  disp(['z0: ', num2str(z0)]);
  disp(['zR: ', num2str(zR)]);

  z_mm = (0:size(bscan1,1)-1)' * dz_mm;
  z0_mm = z0 * dz_mm;
  zR_mm = zR * dz_mm;

  disp(['z0 (mm): ', num2str(z0_mm)]);
  disp(['apparent zR (mm): ', num2str(zR_mm)]);
  disp(['zR (mm): ', num2str(zR_mm/2/1.4)]);

  bscan1_dB = intensity2dB( bscan1 );
  figure; imshownice( bscan1_dB );
  muFit = muFit2D_DRC( bscan2, z_mm, z0_mm, zR_mm, noisePower );
  figure; imshow( muFit, [0 5] );
end

