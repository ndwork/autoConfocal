
function [file1,file2,dz_mm,noisePower] = loadDataCase( datacase )

  mainDir = '/Volumes/ndwork16GB/octData/';
  mainDir = '/Volumes/Seagate2TB/Data/OCTdata/autoConfocal/';

  switch datacase

    case 1
      %dataDir = '/Volumes/Seagate2TB/Data/OCTdata/autoConfocal/160106/phantom_2/';
      dataDir = '/Volumes/ndwork16GB/octData/20160106_focalPlane/phantom_1/';
      file1 = [dataDir, 'phantom_raw_4.raw'];
      file2 = [dataDir, 'phantom_raw_11.raw'];
      dz_mm = 2.57/512;
      noisePower = (1d-2)^2;

    case 2
      dataDir = '/Volumes/ndwork16GB/octData/20160122_verticalOnly/';
      file1 = [dataDir, 'phantom_raw_4.raw'];
      file2 = [dataDir, 'phantom_raw_11.raw'];
      dz_mm = 2.57/512;
      noisePower = (1d-2)^2;

    case 3
      dataDir = '/Volumes/ndwork16GB/octData/20160401_ZeissData/';
      file1 = [dataDir, 'P60698487_HD 5 Line Raster_3-31-2016_13-8-52_OD_sn0057_lineEnhanced.img'];
      file2 = [dataDir, 'P60698487_HD 5 Line Raster_3-31-2016_13-9-10_OD_sn0056_lineEnhanced.img'];
      dz_mm = 2.57/512;
      noisePower = (1d-2)^2;

  end

end
