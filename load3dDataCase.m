
function [vol1_dB,vol2_dB,dz_mm,noisePower] = load3dDataCase( datacase )

  options.dataDimension = 3;
  options.saveInterferogram = 0;

  switch datacase
    case 1
      dataDir = '/Users/ndwork/Dropbox/Sharing/forGenna/autoConfocal/data/';
      datafile1 = [dataDir,'phantom_raw3D_00.raw'];
      datafile2 = [dataDir,'phantom_raw3D_10.raw'];    
      dz_mm = 2.57/512;
      noisePower = (0.5d-2)^2;
  end

  disp('Loading volume 1:');
  [interf1,info1] = getInterferograms(datafile1,options);
  vol1_dB = getBScans(interf1);
  disp('Loading volume 2:');
  [interf2,info2] = getInterferograms(datafile2,options);
  vol2_dB = getBScans(interf2);

end

