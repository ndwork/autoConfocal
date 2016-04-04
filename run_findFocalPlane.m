
function run_findFocalPlane
  clear; close all; rng(1);
  addpath(genpath(pwd))

  %dataDir = '/Volumes/Seagate2TB/Data/OCTdata/autoConfocal/160106/phantom_2/';
  %dataDir = '/Volumes/ndwork16GB/octData/20160106_focalPlane/phantom_1/';
  %dataDir = '/Volumes/ndwork16GB/octData/20160122_verticalOnly/';
  file1 = [dataDir, 'phantom_raw_4.raw'];
  file2 = [dataDir, 'phantom_raw_11.raw'];
  
  [file1,file2,dz_mm,noisePower] = loadDataCase( 3 );

  %[lambda,deltaLambda,dLambda] = getTelestoFalloffParams();

  options.dataDimension = 2;
  options.saveInterferogram = 0;

  [interf1,info1] = getInterferograms(file1,options);
  bscans1_dB = getBScans(interf1);

  [interf2,info2] = getInterferograms(file2,options);
  bscans2_dB = getBScans(interf2);

  tic
  %[z0,zR] = findFocalPlane( bscans1_dB, bscans2_dB );
load junk.mat
  timeTaken = toc;
  disp(['Time taken: ', num2str(timeTaken)]);

  I1 = 10.^( bscans2_dB / 10 );
  z0_mm = z0 * dz_mm;
  zR_mm = zR * dz_mm;
  z_mm = (0:size(I1,2)-1) * dz_mm;
  if exist( 'lambda' ) && exist( 'deltaLambda' )
    muFit = muFit2D_DRC( I1, z_mm, z0_mm, zR_mm, noisePower, ...
      lambda, deltaLambda, dLambda );
  else
    muFit = muFit2D_DRC( I1, z_mm, z0_mm, zR_mm, noisePower );
  end

  figure; imshow( bscans1_dB, [] );  title('Bscan 1');
  figure; imshow( bscans2_dB, [] );  title('Bscan 2');
  figure; imshownice( muFit, 0.5 );

end

