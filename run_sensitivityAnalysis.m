
function run_sensitivityAnalysis
  clear; close all; rng(1);
  addpath(genpath(pwd))

  % Parameters
  mainDir = '/Volumes/ndwork16GB/octData/20170912_sensitivityAnalysis';
  imgIndxs = [3,12];
  ts = -3:3;
  rs = -10:5:10;
  dz_mm = 2.57/512;
  refTransIndx = 4;
  refRotIndx = 1;

  options.dataDimension = 2;
  options.saveInterferogram = 0;

  refDir = [mainDir, '/translation_', num2str(ts(refTransIndx)), ...
    '/angle_', num2str(rs(refRotIndx))];
  files = dir( [refDir,'/OCTData_*.raw'] );
  refFile = [ refDir, '/', files(imgIndxs(1)).name ];
  interf = getInterferograms(refFile,options);
  refBscan_dB = getBScans(interf);
  refBscan = intensity2dB( refBscan_dB, -1 );

  z0s = zeros( numel(ts), numel(rs) );
  zRs = zeros( numel(ts), numel(rs) );

  for tIndx = 1:numel(ts)
   for rIndx = 1:numel(rs)
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

      [z0,zR] = findConfocal_yShearAndTrans( refBscan, thisBscan, ...
        [], [], dz_mm );
      close all;
      z0s(tIndx,rIndx) = z0;
      zRs(tIndx,rIndx) = zR;

    end
  end

end

