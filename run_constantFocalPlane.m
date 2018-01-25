

function run_constantFocalPlane
  clear; close all; rng(1); addpath(genpath(pwd))

  %refIndx = 2;
  %transIndx = 6;
  refIndx = 2;
  transIndx = 6;

  locs = [1 2 3 4 5 7 8 9];
  nLocs = numel( locs );

  trueZ0_mm = 1.5;   % mm
  trueZr_mm = 0.21/2;  % mmOCTData_0.raw
  n = 1.4;  % index of refraction in sample
  trueZR_mm = trueZr_mm * 2 * n;
  Dx = 2.57;  % horiztonal FOV in mm
  Dz = 2.57;  % vertical FOV in mm
  noisePower = (1d5)^2;

  [lambda0, deltaLambda, dLambda] = getTelestoFalloffParams();
  
  mainDir = '/Volumes/ndwork128GB/octData/20171215_responseToReviewers/constantFocalPlane/';
  options.dataDimension = 2;
  options.saveInterferogram = 0;

  hs = cell( nLocs, 1 );
  z0s = zeros( nLocs, 1 );
  zRs = zeros( nLocs, 1 );
  for locIndx = 1:nLocs
  	loc = locs( locIndx );
    thisDir = [ mainDir, 'location', num2str(loc) ];

    refFile = getNthRawFile( [thisDir,'/reference'], refIndx );
    transFile = getNthRawFile( [thisDir,'/transformed'], transIndx );

    [interf1,info1] = getInterferograms(refFile,options);                                   %#ok<ASGLU>
    bscan1_dB = getBScans(interf1);
    bscan1 = intensity2dB( bscan1_dB, -1 );

    [interf2,info2] = getInterferograms(transFile,options);                                 %#ok<ASGLU>
    bscan2_dB = getBScans(interf2);    
    bscan2 = intensity2dB( bscan2_dB, -1 );

    % Remove edges (tend to be erroneous)
    bscan1 = bscan1(1:490,10:end-10);
    bscan2 = bscan2(1:490,10:end-10);
    
    dx_mm = Dx/size(bscan1_dB,2);
    dz_mm = Dz/size(bscan1_dB,1);

    [z0,zR] = findConfocal_yShearAndTrans( bscan1, bscan2, ...
      lambda0, deltaLambda, dLambda, trueZ0_mm, trueZR_mm, dx_mm, dz_mm );

    z_mm = (0:size(bscan1,1)-1)' * dz_mm;
    z0_mm = z0 * dz_mm;  z0s(locIndx) = z0_mm;
    zR_mm = zR * dz_mm;  zRs(locIndx) = zR_mm;

    h = makeConfocalFunction( z_mm, z0_mm, zR_mm );
    hs{locIndx} = h;  
    figure; plotnice( z_mm, h );

    muFit1 = muFit2D_DRC( bscan1, z_mm, z0_mm, zR_mm, noisePower );
    figure; imshowscale( muFit1, 'range', [0 5] );

    close all;
  end

  figure;  hold all;  % Show confocal functions
  for locIndx = 1:nLocs
    plotnice( z_mm, hs{locIndx} );
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




