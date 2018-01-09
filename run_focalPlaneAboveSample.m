
function run_focalPlaneAboveSample
  clear; close all; rng(1); addpath(genpath(pwd))

  %refIndx = 2;
  %transIndx = 6;
  refIndx = 2;
  transIndx = 8;
  nLocs = 10;


  trueZ0_mm = -0.5;   % mm
  trueZr_mm = 0.21/2;  % mmOCTData_0.raw
  n = 1.4;  % index of refraction in sample
  trueZR_mm = trueZr_mm * 2 * n;
  Dx = 2.57;  % horiztonal FOV in mm
  Dz = 2.57;  % vertical FOV in mm
  noisePower = (1d5)^2;

  mainDir = '/Volumes/ndwork128GB/octData/20171215_responseToReviewers/focalPlaneAboveSample/';
  options.dataDimension = 2;
  options.saveInterferogram = 0;

  hs = cell( nLocs, 1 );
  z0s = zeros( nLocs, 1 );
  zRs = zeros( nLocs, 1 );
  for loc = 1:nLocs
    thisDir = [ mainDir, 'location', num2str(loc) ];

    refFile = getNthRawFile( [thisDir,'/reference'], refIndx );
    transFile = getNthRawFile( [thisDir,'/transformed'], transIndx );

    [interf1,info1] = getInterferograms(refFile,options);                                   %#ok<ASGLU>
    bscan1_dB = getBScans(interf1);
    bscan1 = intensity2dB( bscan1_dB, -1 );

    [interf2,info2] = getInterferograms(transFile,options);                                 %#ok<ASGLU>
    bscan2_dB = getBScans(interf2);    
    bscan2 = intensity2dB( bscan2_dB, -1 );

    dx_mm = Dx/size(bscan1_dB,2);
    dz_mm = Dz/size(bscan1_dB,1);

    [z0,zR] = findConfocal_yShearAndTrans( bscan1, bscan2, ...
      trueZ0_mm, trueZR_mm, dx_mm, dz_mm );

    z_mm = (0:size(bscan1,1)-1)' * dz_mm;
    z0_mm = z0 * dz_mm;
    zR_mm = zR * dz_mm;

    h = makeConfocalFunction( z_mm, z0_mm, zR_mm );
    hs{loc} = h;  z0s(loc) = z0_mm;  zRs(loc) = zR_mm;
    figure; plotnice( z_mm, h );

    muFit1 = muFit2D_DRC( bscan1, z_mm, z0_mm, zR_mm, noisePower );
    figure; imshowscale( muFit1, 'range', [0 5] );

    close all;
  end

  figure;  hold all;  % Show confocal functions
  for loc=1:nLocs
    plotnice( z_mm, hs{loc} );
  end

end


function out = getNthRawFile( fileDir, n )
  files = dir( [fileDir,'/*.raw'] );
  nFiles = numel(files);
  
  nums = zeros( nFiles, 1 );
  for i=1:nFiles
    filename = files(i).name;
    nameParts = strsplit( filename, '_' );
    lastParts = strsplit( nameParts{end}, '.' );
    nums(i) = str2double( lastParts{1} );
  end
  [~,sIndxs] = sort( nums );
  files = files( sIndxs );
  out = [ fileDir, '/', files(n).name ];
end


