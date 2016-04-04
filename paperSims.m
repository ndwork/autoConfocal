
function paperSims
  clear; close all; rng(1);
  addpath(genpath(pwd))

  outDir = '.';
  %mkdir(outDir);

  fid = fopen([outDir,'/paperSimsOut.csv'],'w');
  fprintf(fid, ['name, algorithm, z0, z0_true, zR, zR_true, mus, thicks, ', ...
    'noiseProportion, meanEtbDepth, time taken (s) \n']);

  noiseProportion = 1d-5;
  noisePower = (0.5d-2)^2;
  N = 100;
  simIndx = 0;



  %%  Focal Length Error
  disp('Analyzing Focal length errors');
  fprintf( fid, 'Focal length error analysis  \n');

  z0_true = 0.5;

  simMus{1} = [1,2];
  simThicks{1} = [1];

  zR = 0.1059 * 2 * 1.37;   % n = 1.37

  for j=1:numel(simMus)

    muStrings = strtrim(cellstr(num2str(simMus{j}'))');
    muString = strjoin( muStrings, '_' );

    [I, z, z0, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
      N, -1, z0_true, noiseProportion, simMus{j}, [1], zR );
    dz = z(2) - z(1);

    thickStrings = strtrim(cellstr(num2str(simThicks{1}'))');
    thickStrings{end+1} = num2str( max(z) - sum(simThicks{1}) );
    thickString = strjoin( thickStrings, '_' );

    for z0=z0_true-0.1:0.01:z0_true+0.1;
      simIndx = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
        z0_true, zR, muString, thickString, noiseProportion, noisePower, trueMu );
    end
  end
  clear simMus;
  clear simThicks;



  %%  Rayleigh Range Error
  disp('Analyzing Rayleigh range errors');
  fprintf( fid, 'Rayleigh range error analysis  \n');

  z0 = 0.5;

  simMus{1} = [1,2];
  simThicks{1} = [1];

  zR_true = 0.1059 * 2 * 1.37;   % n = 1.37

  for j=1:numel(simMus)

    muStrings = strtrim(cellstr(num2str(simMus{j}'))');
    muString = strjoin( muStrings, '_' );

    [I, z, z0, zR, muAlpha, muBeta, muL0, trueMu] = makePhantom2D( ...
      N, -1, z0, noiseProportion, simMus{j}, [1], zR_true );
    dz = z(2) - z(1);

    thickStrings = strtrim(cellstr(num2str(simThicks{1}'))');
    thickStrings{end+1} = num2str( max(z) - sum(simThicks{1}) );
    thickString = strjoin( thickStrings, '_' );

    %for zR=0:zR_true/10:2*zR_true;
    for zR = zR_true-0.1 : 0.1/20 : zR_true+0.1
      simIndx = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
        z0, zR_true, muString, thickString, noiseProportion, noisePower, trueMu );
    end
  end
  clear simMus;
  clear simThicks;



  %% Cleanup
  fclose(fid);
end



function simIndx_out = evalCodes( simIndx, I, z, z0, zR, fid, outDir, ...
  z0_true, zR_true, muString, thickString, noiseProportion, noisePower, trueMu )

  offsetThreshPercent = 0.05;

  tic;
  muFit = muFit2D_DRC( I, z, z0, zR, noisePower );
  timeTaken = toc;
  meanEtbDepth = findErrorMetrics( muFit, trueMu, z, offsetThreshPercent );
  filename = ['sim_', num2str(simIndx,'%5.5i'),'.mat'];
  fprintf( fid, [filename, ', DRC, %8.3f, %8.3f, %8.3f, %8.3f, ', muString, ...
        ', ', thickString, ', %5.2e, %8.3f, %16.11f\n'], ...
        z0, z0_true, zR, zR_true, noiseProportion, meanEtbDepth, timeTaken );
  simIndx = simIndx + 1;
  
  simIndx_out = simIndx;
end


