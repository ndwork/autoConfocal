
function [bscan1,bscan2,dz_mm,noisePower,trans] = loadDataCase( datacase )

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
      dataDir = [mainDir,'20160106_focalPlane/phantom_1/'];
      file1 = [dataDir, 'phantom_raw_4.raw'];
      file2 = [dataDir, 'phantom_raw_11.raw'];
      dz_mm = 2.57/512;
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      noisePower = (0.5d-2)^2;
      trans = 'vShift';

    case 2
      dataDir = [mainDir,'/20160122_verticalOnly/'];
      file1 = [dataDir, 'phantom_4.raw'];
      %file2 = [dataDir, 'phantom_11.raw'];
      file2 = [dataDir, 'phantom_18.raw'];
      dz_mm = 2.57/512;
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      noisePower = (0.5d-2)^2;
      trans = 'vShift';

    case 3
      dataDir = [mainDir,'/20160401_ZeissData/'];
      file1 = [dataDir, 'P60698487_HD 5 Line Raster_3-31-2016_13-9-32_OD_sn0058_lineEnhanced.img'];
      file2 = [dataDir, 'P60698487_HD 5 Line Raster_3-31-2016_13-9-10_OD_sn0056_lineEnhanced.img'];
      %file1 = [dataDir, 'P60698487_HD 5 Line Raster_3-31-2016_13-8-52_OD_sn0057_lineEnhanced.img'];
      dz_mm = 2.00/1024;
      bscan1_dB = loadLengBscan(file1);
      bscan1_dB = bscan1_dB(:,30:end,1);
      bscan1_dB = rot90(bscan1_dB,1);
      bscan2_dB = loadLengBscan(file2);
      bscan2_dB = bscan2_dB(:,30:end,1);
      bscan2_dB = rot90(bscan2_dB,1);
      noisePower = (0.5d-2)^2;
      trans = 'yShearAndTrans';

    case 4
      dataDir = [mainDir,'/20160418_angles/'];
      file1 = [dataDir, 'negativeAngle/16020404_17.raw'];
      file2 = [dataDir, 'flat/16020404_31.raw'];
      dz_mm = 2.57/512;
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      noisePower = (0.5d-2)^2;

    case 5
      % structured phantom with rotation and translation
      dataDir = [mainDir,'/20160818_layeredPhantom/'];
      file1 = [dataDir,'phantom_raw_0.raw'];
      file2 = [dataDir,'phantom_raw_9.raw'];
      dz_mm = 2.57/512;
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      noisePower = (0.5d-2)^2;
      trans = 'yShearAndTrans';

    case 6
      % structured phantom with vertical shift
      dataDir = [mainDir,'/20160818_layeredPhantom/'];
      file1 = [dataDir,'phantom_raw_0.raw'];
      file2 = [dataDir,'phantom_raw_3.raw'];
      dz_mm = 2.57/512;
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      noisePower = (0.5d-2)^2;
      trans = 'vShift';

    case 7
      % structured phantom with vertical shift
      dataDir = [mainDir,'/20160818_layeredPhantom/'];
      file1 = [dataDir,'phantom_raw_14.raw'];
      file2 = [dataDir,'phantom_raw_18.raw'];
      dz_mm = 2.57/512;
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      noisePower = (0.5d-2)^2;
      trans = 'vShift';
      
    case 8
      %structured phantom with vertical shift
      dataDir = [mainDir,'/20160808_uniformPhantom/'];
      file1 = [dataDir,'phantom_raw_20.raw'];
      file2 = [dataDir,'phantom_raw_26.raw'];
      dz_mm = 2.57/512;
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      noisePower = (0.5d-2)^2;
      trans = 'vShift';

    case 9
      % structured phantom with vertical shift
      dataDir = ['/Users/ndwork/Dropbox/Sharing/forGenna/20160922_layeredPhantom/set3/'];
      file1 = [dataDir,'phantom_raw_27.raw'];
      file2 = [dataDir,'phantom_raw_29.raw'];
      dz_mm = 2.57/512;
      [interf1,info1] = getInterferograms(file1,options);
      bscan1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscan2_dB = getBScans(interf2);
      noisePower = (0.5d-2)^2;
      trans = 'vShift';

  end

  bscan1 = intensity2dB( bscan1_dB, -1 );
  bscan2 = intensity2dB( bscan2_dB, -1 );
end
