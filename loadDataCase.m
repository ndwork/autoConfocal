
function [bscans1_dB,bscans2_dB,dz_mm,noisePower] = loadDataCase( datacase )

  mainDir = '/Volumes/ndwork16GB/octData/';
  %mainDir = '/Volumes/Seagate2TB/Data/OCTdata/autoConfocal/';

  options.dataDimension = 2;
  options.saveInterferogram = 0;

  switch datacase

    case 1
      dataDir = [mainDir,'20160106_focalPlane/phantom_1/'];
      file1 = [dataDir, 'phantom_raw_4.raw'];
      file2 = [dataDir, 'phantom_raw_11.raw'];
      dz_mm = 2.57/512;
      noisePower = (1d-2)^2;
      [interf1,info1] = getInterferograms(file1,options);
      bscans1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscans2_dB = getBScans(interf2);

    case 2
      dataDir = [mainDir,'/20160122_verticalOnly/']';
      file1 = [dataDir, 'phantom_raw_4.raw'];
      file2 = [dataDir, 'phantom_raw_11.raw'];
      dz_mm = 2.57/512;
      noisePower = (1d-2)^2;
      [interf1,info1] = getInterferograms(file1,options);
      bscans1_dB = getBScans(interf1);
      [interf2,info2] = getInterferograms(file2,options);
      bscans2_dB = getBScans(interf2);
      
      

    case 3
      dataDir = [mainDir,'/20160401_ZeissData/'];
      file1 = [dataDir, 'P60698487_HD 5 Line Raster_3-31-2016_13-8-52_OD_sn0057_lineEnhanced.img'];
      file2 = [dataDir, 'P60698487_HD 5 Line Raster_3-31-2016_13-9-10_OD_sn0056_lineEnhanced.img'];
      dz_mm = 2.57/512;
      noisePower = (1d-2)^2;
      bscans1_dB = loadLengBscan(file1);
      bscans1_dB = bscans1_dB(:,:,1);
      bscans1_dB = rot90(bscans1_dB,1);
      bscans2_dB = loadLengBscan(file2);
      bscans2_dB = bscans2_dB(:,:,1);
      bscans2_dB = rot90(bscans2_dB,1);

  end

end
