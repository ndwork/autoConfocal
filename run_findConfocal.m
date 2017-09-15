
function run_findConfocal
  clear; close all; rng(1);
  addpath(genpath(pwd))

  %datacase =  8;  % Structured phantom with rotation and translation
  datacase = 11;  % Zeiss retina
  %datacase = 12;  % Rabbit eye
  %datacase = 14;  % Sensitivity analysis phantom
  [bscan1,bscan2,dz_mm,noisePower,trans,trueZ0_mm,trueZR_mm] = ...
    loadDataCase( datacase );

  switch trans
    case 'vShift'
      [z0,zR] = findConfocal_vShift( bscan1, bscan2, ...
        trueZ0_mm, trueZR_mm, dz_mm );
    case 'vShiftRot'
      [z0,zR] = findConfocal_vShiftRot( bscan1, bscan2 );
    case 'yShearAndTrans'
      [z0,zR] = findConfocal_yShearAndTrans( bscan1, bscan2, ...
        trueZ0_mm, trueZR_mm, dz_mm );
  end

  disp(['z0: ', num2str(z0)]);
  disp(['zR: ', num2str(zR)]);

  z_mm = (0:size(bscan1,1)-1)' * dz_mm;
  z0_mm = z0 * dz_mm;
  zR_mm = zR * dz_mm;

  disp(['z0 (mm): ', num2str(z0_mm)]);
  disp(['zR (mm): ', num2str(zR_mm)]);
  disp(['True z0 (mm): ', num2str(trueZ0_mm) ]);
  disp(['True zR (mm): ', num2str(trueZR_mm) ]);

  bscan1_dB = intensity2dB( bscan1 );
  figure; imshownice( bscan1_dB );
  title('B-Scan 1 in dB');


  tic;
  muFit = muFit2D_DRC( bscan1, z_mm, z0_mm, zR_mm, noisePower );
  toc
  figure; imshow( muFit, [0 2.5] );

  if datacase == 11
    % find the attenuation coefficient of the nerve fiber layer
    mask = muFit > 2;
    mask = imdilate( mask, strel( 'square', 8 ) );
    mask = imerode( mask, strel( 'square', 10 ) );
    regions = bwlabel( mask );
    nflRegionIndx = regions(348,20);
    nflMus = muFit( regions == nflRegionIndx );
    disp(['NFL Mean Mu: ', num2str(mean(nflMus))]);
    disp(['NFL Median Mu: ', num2str(median(nflMus))]);
  end

  if datacase == 14
    % find the median attenuation coefficient of the first and second
    % layers
    segLayer1 = muFit( 67:117, 18:224 );
    med1 = median(segLayer1(:));
    disp(['Median of layer 1: ', num2str(med1)]);
    segLayer2 = muFit(137:225,266:384);
    med2 = median(segLayer2(:));
    disp(['Median of layer 2: ', num2str(med2)]);
  end
end

