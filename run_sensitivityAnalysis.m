
function run_sensitivityAnalysis
  clear; close all; rng(1);
  addpath(genpath(pwd))

  % Parameters
  mainDir = '/Volumes/ndwork16GB/octData/20170912_sensitivityAnalysis';
  imgIndxs = [3,12];
  ts = -3:3;
  rs = -10:5:10;
  dx_mm = 2.57/512;
  dz_mm = 2.57/512;
  refTransIndx = 4;
  refRotIndx = 1;
  trueZ0_mm = 1.5;
  trueZR_mm = 0.28;

  options.dataDimension = 2;
  options.saveInterferogram = 0;

  refDir = [mainDir, '/translation_', num2str(ts(refTransIndx)), ...
    '/angle_', num2str(rs(refRotIndx))];
  files = dir( [refDir,'/OCTData_*.raw'] );
  refFile = [ refDir, '/', files(imgIndxs(1)).name ];
  interf = getInterferograms(refFile,options);
  refBscan_dB = getBScans(interf);
  refBscan = intensity2dB( refBscan_dB, -1 );

  z0Cells = cell( 1, numel(ts) );
  zRCells = cell( 1, numel(ts) );
  overlapCells = cell( 1, numel(ts) );
  parfor tIndx = 1:numel(ts)
    theseZ0s = zeros(numel(rs),1);
    theseZRs = zeros(numel(rs),1);
    theseOverlaps = zeros(numel(rs),1);

    for rIndx = 1:numel(rs)
      rng(2);
      disp([ 'Working on Translation ', num2str(tIndx), ' of ', num2str(numel(ts)) ]);
      disp([ 'Working on Rotation ', num2str(rIndx), ' of ', num2str(numel(rs)) ]);

      thisTrans = ts(tIndx);
      thisRot = rs(rIndx);

      thisDir = [mainDir, '/translation_', num2str(thisTrans), ...
        '/angle_', num2str(thisRot)];
      files = dir( [thisDir,'/OCTData_*.raw'] );
      thisFile = [ thisDir, '/', files(imgIndxs(2)).name ];

      interf = getInterferograms(thisFile,options);
      thisBscan_dB = getBScans(interf);
      thisBscan = intensity2dB( thisBscan_dB, -1 );

      [z0,zR,overlap] = findConfocal_yShearAndTrans( refBscan, thisBscan, ...
        [], [], dx_mm, dz_mm );
      close all;
      theseZ0s(rIndx) = z0;
      theseZRs(rIndx) = zR;
      theseOverlaps(rIndx) = overlap;
    end

    z0Cells{tIndx} = theseZ0s;
    zRCells{tIndx} = theseZRs;
    overlapCells{tIndx} = theseOverlaps;
  end

  z0s = cell2mat( z0Cells );  z0s_mm = z0s * dz_mm;
  zRs = cell2mat( zRCells );  zRs_mm = zRs * dz_mm;
  overlaps_mmSq = cell2mat( overlapCells );

  disp('z0s (mm):');  disp(z0s_mm);
  disp('zRs (mm):');  disp(zRs_mm);
  disp('Overlaps (mm^2):');  disp(overlaps_mmSq);

  errorsZ0_mm = trueZ0_mm - z0s_mm;
  errorsZR_mm = trueZR_mm - zRs_mm;

  [overlaps_mmSq,sortedIndxs] = sort(overlaps_mmSq(:));
  errorsZ0_mm = errorsZ0_mm(sortedIndxs);
  errorsZR_mm = errorsZR_mm(sortedIndxs);
  figure; plotnice( overlaps_mmSq(:), abs(errorsZ0_mm) );  title('errorsZ0_mm');
  figure; plotnice( overlaps_mmSq(:), abs(errorsZR_mm) );  title('errorsZR_mm');
end

