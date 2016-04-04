
function run_findFocalPlane
  clear; close all; rng(1);
  addpath(genpath(pwd))

  datacase = 3;

  [bscans1_dB,bscans2_dB,dz_mm,noisePower] = loadDataCase( datacase );

  %[lambda,deltaLambda,dLambda] = getTelestoFalloffParams();


  tic
  [z0,zR] = findFocalPlane( bscans1_dB, bscans2_dB );
%load junk.mat
  timeTaken = toc;
  disp(['Time taken: ', num2str(timeTaken)]);

  I1 = 10.^( bscans2_dB / 10 );
  z0_mm = z0 * dz_mm;
  zR_mm = zR * dz_mm;
  z_mm = (0:size(I1,2)-1) * dz_mm;
  if exist( 'lambda', 'var' ) && exist( 'deltaLambda', 'var' )
    muFit = muFit2D_DRC( I1, z_mm, z0_mm, zR_mm, noisePower, ...
      lambda, deltaLambda, dLambda );
  else
    muFit = muFit2D_DRC( I1, z_mm, z0_mm, zR_mm, noisePower );
  end

  figure; imshow( bscans1_dB, [] );  title('Bscan 1');
  figure; imshow( bscans2_dB, [] );  title('Bscan 2');
  figure; imshownice( muFit, 0.5 );

end

