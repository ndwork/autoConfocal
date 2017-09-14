
function [bscan1,bscan2,dz_mm,noisePower,trans,trueZ0_mm,trueZR_mm] = ...
  loadDataCase( datacase )

  if exist( '/Volumes/Seagate2TB/Data/OCTdata/autoConfocal/', 'dir' )
    mainDir = '/Volumes/Seagate2TB/Data/OCTdata/autoConfocal/';
  elseif exist( '/Volumes/ndwork16GB/octData/', 'dir' ) 
    mainDir = '/Volumes/ndwork16GB/octData/';
  else
    error('Could not find main directory');
  end

  options.dataDimension = 2;
  options.saveInterferogram = 0;

  switch datacase

    case 1
      % structured phantom with vertical shift
      trueZ0_mm = 0.86;  % mm
      trueZR_mm = 0.29;  % mm
      file1 = [mainDir,'/20160922_layeredPhantom/set3_withoutHoles/phantom_raw_27.raw'];
      file2 = [mainDir,'/20160922_layeredPhantom/set3_withoutHoles/phantom_raw_29.raw'];
      dz_mm = 2.57/512;
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      noisePower = (0.5d-2)^2;
      trans = 'vShift';

    case 2
      % structured phantom with vertical shift
      trueZ0_mm = 1.24;  % mm
      trueZR_mm = 0.29;  % mm
      file1 = [mainDir,'/20160922_layeredPhantom/set4/flat/phantom_raw_0.raw'];
      file2 = [mainDir,'/20160922_layeredPhantom/set4/flat/phantom_raw_5.raw'];
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      dz_mm = 2.57/size(bscan1_dB,1);
      noisePower = (0.5d-2)^2;
      trans = 'vShift';

    case 3
      % structured phantom with vertical shift
      trueZ0_mm = 1.35;   % mm
      trueZR_mm = 0.29;  % mm
      file1 = [mainDir,'/20160922_layeredPhantom/set5/phantom_raw_41.raw'];
      file2 = [mainDir,'/20160922_layeredPhantom/set5/phantom_raw_45.raw'];
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      dz_mm = 2.57/size(bscan1_dB,1);
      noisePower = (0.5d-2)^2;
      trans = 'vShift';

    case 4
      % structured phantom with shift and rotation
      trueZ0_mm = 1.7;   % mm
      trueZR_mm = 0.29;  % mm
      file1 = [mainDir,'/20160922_layeredPhantom/set5/phantom_raw_41.raw'];
      file2 = [mainDir,'/20160922_layeredPhantom/set5/phantom_raw_34.raw'];
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      dz_mm = 2.57/size(bscan1_dB,1);
      %noisePower = (0.5d-2)^2;
      %noisePower = (1d-1)^2;
      noisePower = 10000000;
      trans = 'yShearAndTrans';

    case 5
      % structured phantom with shift and rotation
      file1 = [mainDir,'/20161003_layeredPhantom/set1/phantom_raw_16.raw'];
      file2 = [mainDir,'/20161003_layeredPhantom/set1/phantom_raw_07.raw'];
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      dz_mm = 2.57/size(bscan1_dB,1);
      noisePower = (0.5d-2)^2;
      trans = 'yShearAndTrans';
      
    case 6
      % structured phantom with shift and rotation
      trueZ0_mm = 1.35;   % mm
      trueZR_mm = 0.29;  % mm
      file1 = [mainDir,'/20161003_layeredPhantom/set2/phantom_raw_01.raw'];
      file2 = [mainDir,'/20161003_layeredPhantom/set2/phantom_raw_21.raw'];
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      dz_mm = 2.57/size(bscan1_dB,1);
      noisePower = (0.5d-2)^2;
      %noisePower = (1d-1)^2;
      trans = 'yShearAndTrans';

    case 7
      % structured phantom with shift and rotation
      trueZ0_mm = 1.03;  % mm
      trueZR_mm = 0.29;  % mm
      file1 = [mainDir,'/20161003_layeredPhantom/set3/phantom_raw_00.raw'];
      file2 = [mainDir,'/20161003_layeredPhantom/set3/phantom_raw_16.raw'];
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      %bscan1_dB = bscan1_dB(1:400,:);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      dz_mm = 2.57/size(bscan1_dB,1);
      %bscan2_dB = bscan2_dB(1:400,:);
      %noisePower = (0.5d-2)^2;
      noisePower = (1d-1)^2;
      trans = 'yShearAndTrans';

    case 8
      % structured phantom with shift and rotation
      trueZ0_mm = 0.9;   % mm
      trueZR_mm = 0.29;  % mm
      file1 = [mainDir,'/20161005_layeredPhantom/set1/phantom_raw_1.raw'];
      file2 = [mainDir,'/20161005_layeredPhantom/set1/phantom_raw_20.raw'];
      %file2 = [mainDir,'/20161005_layeredPhantom/set1/phantom_raw_23.raw'];
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      %bscan1_dB = bscan1_dB(1:400,:);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      dz_mm = 2.57/size(bscan1_dB,1);
bscan1_dB = bscan1_dB(1:400,:);
bscan2_dB = bscan2_dB(1:400,:);
      %bscan2_dB = bscan2_dB(1:400,:);
      %noisePower = (0.5d-2)^2;
      %noisePower = 100000000;
      noisePower = 1;
      trans = 'yShearAndTrans';

    case 9
      % Zeiss data from Dr. Leng
      % eye with shift and rotation
      trueZ0_mm = 0;   % mm
      trueZR_mm = 0;  % mm
      file1 = [mainDir,'/20160401_ZeissData/P60698487_HD 5 Line Raster_3-31-2016_13-9-32_OD_sn0058_lineEnhanced.img'];
      file2 = [mainDir,'/20160401_ZeissData/P60698487_HD 5 Line Raster_3-31-2016_13-9-10_OD_sn0056_lineEnhanced.img'];
      %file2 = [mainDir,'/20160401_ZeissData/P60698487_HD 5 Line Raster_3-31-2016_13-8-52_OD_sn0057_lineEnhanced.img'];
      tmp1 = loadLengBscan(file1);
      bscan1_dB = rot90( tmp1(:,:,2) );
      tmp2 = loadLengBscan(file2);
      bscan2_dB = rot90( tmp2(:,:,2) );
      dz_mm = 2.00/size(bscan1_dB,1);

      % determine scale and shift parameters
      [k,c] = findZeissScaleShift( bscan1_dB, dz_mm );
      bscan1_dB = 10.^( (bscan1_dB-c) / k );
      bscan2_dB = 10.^( (bscan2_dB-c) / k );

      noisePower = 10000000;
      trans = 'yShearAndTrans';

    case 10
      % Zeiss data from Dr. Leng
      % structured phantom with shift and rotation
      trueZ0_mm = 0;   % mm
      trueZR_mm = 0;  % mm
      file1 = [mainDir,'/20160421_ZeissData/phantom/PCZMI371589661_HD 5 Line Raster_4-21-2016_13-16-48_OD_sn0078_lineEnhanced.img'];
      file2 = [mainDir,'/20160421_ZeissData/phantom/PCZMI371589661_HD 5 Line Raster_4-21-2016_13-17-3_OD_sn0079_lineEnhanced.img'];
      %file2 = [mainDir,'/20160421_ZeissData/phantom/PCZMI371589661_HD 5 Line Raster_4-21-2016_13-17-31_OD_sn0080_lineEnhanced.img'];
      tmp1 = loadLengBscan(file1);
      bscan1_dB = rot90( tmp1(:,:,2), 3 );
      tmp2 = loadLengBscan(file2);
      bscan2_dB = rot90( tmp2(:,:,2), 3 );
      dz_mm = 2.00/size(bscan1_dB,1);

      bscan1_dB = bscan1_dB / 4;
      bscan2_dB = bscan2_dB / 4;

      noisePower = (0.5d-2)^2;
      trans = 'vShift';

    case 11
      % Zeiss data from Dr. Leng
      % structured phantom with shift and rotation
      trueZ0_mm = 0;   % mm
      trueZR_mm = 0;  % mm
      file2 = [mainDir,'/20160421_ZeissData/eye/P60698487_HD 5 Line Raster_4-21-2016_13-24-11_OD_sn0084_lineEnhanced.img'];
      file1 = [mainDir,'/20160421_ZeissData/eye/P60698487_HD 5 Line Raster_4-21-2016_13-24-39_OD_sn0085_lineEnhanced.img'];
      %file2 = [mainDir,'/20160421_ZeissData/eye/P60698487_HD 5 Line Raster_4-21-2016_13-24-56_OD_sn0086_lineEnhanced.img'];
      tmp1 = loadLengBscan(file1);
      bscan1_dB = rot90( tmp1(:,:,2) );
      dz_mm = 2.00/size(bscan1_dB,1);
      cropAmount = size(bscan1_dB)-30;
      bscan1_dB = cropData( bscan1_dB, cropAmount );
      tmp2 = loadLengBscan(file2);
      bscan2_dB = rot90( tmp2(:,:,2) );
      bscan2_dB = cropData( bscan2_dB, cropAmount );

      bscan1_dB = bscan1_dB / 4;
      bscan2_dB = bscan2_dB / 4;

      bscan1_dB = bscan1_dB(:,100:end-100);
      bscan2_dB = bscan2_dB(:,100:end-100);

      noisePower = (0.5d-2)^2;
      trans = 'yShearAndTrans';

    case 12
      % Rabbit outer eye
      trueZ0_mm = 0.7;   % mm
      trueZR_mm = 0.29;  % mm
      file2 = [mainDir,'/20170215_rabbit/outerEye/flat/phantom_raw_68.raw'];
      file1 = [mainDir,'/20170215_rabbit/outerEye/rotated/phantom_raw_86.raw'];
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      dz_mm = 2.57/size(bscan1_dB,1);
      %bscan1_dB = bscan1_dB(1:400,:);
      %bscan2_dB = bscan2_dB(1:400,:);
      noisePower = (0.5d-2)^2;
      %noisePower = 100000000;
      trans = 'yShearAndTrans';

    case 13
      % layered phantom with vertical shift
      trueZ0_mm = 1.1;  % mm
      trueZR_mm = 0.29;  % mm
      file1 = [mainDir,'/20160106_focalPlane/phantom_2/phantom_raw_4.raw'];
      file2 = [mainDir,'/20160106_focalPlane/phantom_2/phantom_raw_9.raw'];
      dz_mm = 2.57/512;
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      bscan1_dB = bscan1_dB(1:400,:);
      bscan2_dB = bscan2_dB(1:400,:);
      noisePower = (0.5d-2)^2;
      trans = 'vShift';

    case 14
      % structured phantom with shift and rotation
      trueZ0_mm = 1.2;   % mm
      trueZR_mm = 0.28;  % mm
      file1 = [mainDir,'/20170912_sensitivityAnalysis/translation_0/angle_0/OCTData_03.raw'];
      file2 = [mainDir,'/20170912_sensitivityAnalysis/translation_0/angle_0/OCTData_09.raw'];
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      %bscan1_dB = bscan1_dB(1:400,:);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      dz_mm = 2.57/size(bscan1_dB,1);
%bscan1_dB = bscan1_dB(1:400,:);
%bscan2_dB = bscan2_dB(1:400,:);
      noisePower = (0.5d-2)^2;
      trans = 'yShearAndTrans';
  end

  bscan1 = intensity2dB( bscan1_dB, -1 );
  bscan2 = intensity2dB( bscan2_dB, -1 );
end

