
function run_findConfocal
  clear; close all; rng(1);
  addpath(genpath(pwd))

  %datacase =  1;  % Structured phantom with rotation and translation
  datacase = 11;  % Zeiss retina
  %datacase = 12;  % Rabbit eye
  %datacase = 14;   % Structured phantom 2 (of sensitivity analysis ) 
                   % with rotation and translation
  %datacase = 15;  % Sensitivity analysis phantom
  %datacase = 16;  % Just noise
  [bscan1,bscan2,dx_mm,dz_mm,noisePower,trans,trueZ0_mm,trueZR_mm] = ...
    loadDataCase( datacase );

  switch trans
    case 'vShift'
      [z0,zR] = findConfocal_vShift( bscan1, bscan2, ...
        trueZ0_mm, trueZR_mm, dx_mm, dz_mm );
    case 'vShiftRot'
      [z0,zR] = findConfocal_vShiftRot( bscan1, bscan2 );
    case 'yShearAndTrans'
      [z0,zR] = findConfocal_yShearAndTrans( bscan1, bscan2, ...
        trueZ0_mm, trueZR_mm, dx_mm, dz_mm );
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
  figure; imshownice( bscan1_dB ); title('B-Scan 1 in dB');

  bscan2_dB = intensity2dB( bscan2 );
  figure; imshownice( bscan2_dB ); title('B-Scan 2 in dB');


  tic;
  muFit1 = muFit2D_DRC( bscan1, z_mm, z0_mm, zR_mm, noisePower );
  toc
  figure; imshowscale( muFit1, 3, 'range', [0 5] );
  colormap('jet'); colorbarnice;
  title('muFit of 1st B-scan');

  figure; imshowscale( 20*log10(muFit1), 3, 'range', [-15 15] );
  colormap('jet'); colorbarnice;
  title('muFit of 1st B-scan in dB');

  tic;
  muFit2 = muFit2D_DRC( bscan2, z_mm, z0_mm, zR_mm, noisePower ); colorbarnice;
  toc
  figure; imshowscale( muFit2, 3, 'range', [0 5] );  colormap('jet'); colorbarnice;
  title('mufit of 2nd B-scan');

  figure; imshowscale( 20*log10(muFit2), 3, 'range', [-15 15] );
  colormap('jet'); colorbarnice;
  title('muFit2 of 1st B-scan in dB');

  if datacase == 11
    % find the attenuation coefficient of the nerve fiber layer
    mask = muFit1 > 2;
    mask = imdilate( mask, strel( 'square', 8 ) );
    mask = imerode( mask, strel( 'square', 10 ) );
    regions = bwlabel( mask );
    nflRegionIndx = regions(348,20);
    nflMus = muFit1( regions == nflRegionIndx );
    disp([ 'NFL Mean Mu: ', num2str(mean(nflMus)) ]);
    disp([ 'NFL Median Mu: ', num2str(median(nflMus)) ]);
    disp([ 'NFL Variance of Mus: ', num2str(var(nflMus)) ]);
    disp([ 'NFL Standard Deviation of Mus: ', num2str(std(nflMus)) ]);
  end

  if datacase == 12
    lowY1 = 154;   highY1 = 171;
    lowX1 = 177;   highX1 = 222;

    figure;  imshownice( bscan1_dB, 3 );  title('bscan1_dB');
    features1 = [ [lowX1;lowX1;highX1;highX1], ...
                  [lowY1;highY1;lowY1;highY1]  ];
    showFeaturesOnImg( features1, 'scale', 3, 'color', 'k' );
    
    subMu1 = muFit1( lowY1:highY1, lowX1:highX1 );
    subBScan1 = bscan1( lowY1:highY1, lowX1:highX1 );
    z1_mm = z_mm(lowY1:highY1);
    meanMu1 = mean( subMu1(:) );
    disp([ 'Mean of sclera region 1 mu is: ', num2str(meanMu1) ]);
    stdMu1 = std( subMu1(:) );
    disp([ 'Std Dev of sclera region 1 mu is: ', num2str(stdMu1) ]);

    meanData1 = mean( subBScan1, 2 );
    stdData1 = std( subMu1, [], 2 );
    muFaber = muFitFaber( meanData1, z1_mm, z0_mm, zR_mm, stdData1 );
    disp([ 'Mu of region 1 with Faber Method: ', num2str(muFaber) ]);

    lowY2 = 93;  highY2 = 116;
    lowX2 = 291;  highX2 = 345;
    

    figure;  imshownice( bscan2_dB, 3 );  title('bscan2_dB');
    features2 = [ [lowX2;lowX2;highX2;highX2], ...
                  [lowY2;highY2;lowY2;highY2]  ];
    showFeaturesOnImg( features2, 'scale', 3, 'color', 'k' );

    subMu2 = muFit2( lowY2:highY2, lowX2:highX2 );
    subBScan2 = bscan2( lowY2:highY2, lowX2:highX2 );
    z2_mm = z_mm( lowY2:highY2 );
    meanMu2 = mean( subMu2(:) );
    disp(['Mean of sclera region 2 mu is: ', num2str(meanMu2)]);
    stdMu2 = std( subMu2(:) );
    disp(['Std Dev of sclera region 2 mu is: ', num2str(stdMu2)]);

    meanData2 = mean( subBScan2, 2 );
    stdData2 = std( subMu2, [], 2 );
    muFaber2 = muFitFaber( meanData2, z2_mm, z0_mm, zR_mm, stdData2 );
    disp([ 'Mu of region 2 with Faber Method: ', num2str(muFaber2) ]);
  end

  if datacase == 14
    % find the median attenuation coefficient of the first and second
    % layers
    segLayer1a = muFit1( 67:117, 18:224 );
    segLayer1b = muFit1( 173:263, 430:500 );
    segLayer1 = [ segLayer1a(:); segLayer1b(:); ];
    disp(['Median of layer 1: ', num2str(median(segLayer1(:)))]);
    disp(['Mean of layer 1: ', num2str(mean(segLayer1(:)))]);
    segLayer2 = muFit1(137:225,266:384);
    disp(['Median of layer 2: ', num2str(median(segLayer2(:)))]);
    disp(['Mean of layer 2: ', num2str(mean(segLayer2(:)))]);
    rotMuFit = imrotate( muFit1, 15, 'nearest' );
    segLayer3 = rotMuFit( 322:360, 114:506 );
    disp(['Median of layer 3: ', num2str(median(segLayer3(:)))]);
    disp(['Mean of layer 3: ', num2str(mean(segLayer3(:)))]);
    segLayer4 = rotMuFit( 412:460, 128:490 );
    disp(['Median of layer 4: ', num2str(median(segLayer4(:)))]);
    disp(['Mean of layer 4: ', num2str(mean(segLayer4(:)))]);
    segLayer5 = muFit1( 417:463, 96:208 );
    disp(['Median of layer 5: ', num2str(median(segLayer5(:)))]);
    disp(['Mean of layer 5: ', num2str(mean(segLayer5(:)))]);
  end
end

