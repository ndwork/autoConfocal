
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
      %noisePower = (0.5d-2)^2;
      noisePower = (1d-1)^2;
      trans = 'yShearAndTrans';

    case 7
      % structured phantom with shift and rotation
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
      file2 = [mainDir,'/20161005_layeredPhantom/set1/phantom_raw_23.raw'];
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      %bscan1_dB = bscan1_dB(1:400,:);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      dz_mm = 2.57/size(bscan1_dB,1);
      %bscan2_dB = bscan2_dB(1:400,:);
      %noisePower = (0.5d-2)^2;
      noisePower = 100000000;
      trans = 'yShearAndTrans';
      
  end

  bscan1 = intensity2dB( bscan1_dB, -1 );
  bscan2 = intensity2dB( bscan2_dB, -1 );

  %bscan1 = bscan1 / 55743;
  %bscan2 = bscan2 / 55743;
end
