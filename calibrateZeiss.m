
function calibrateZeiss
  close all;  clear;

  datacase = 10;
  [bscan1,bscan2,dz_mm] = loadDataCase( datacase );
  zeissData1 = intensity2dB( bscan1 );
  zeissData2 = intensity2dB( bscan2 );
  mu = 8;

  for imgIndx=1:2

    if imgIndx==1
      shearAngle = 20;
      colStart = 450;   colEnd = 550;
      rowStart = 280;   rowEnd = 380;
      bscan = bscan1;
      zeissData = zeissData1;
    elseif imgIndx==2
      shearAngle = 28;
      colStart = 470;   colEnd = 530;
      rowStart = 490;   rowEnd = 610;
      bscan = bscan2;
      zeissData = zeissData2;
    end

    shearedImg = shearImg( zeissData, shearAngle*pi/180 );
    dataValues = mean( shearedImg(rowStart:rowEnd,colStart:colEnd), 2 );
    zs = (rowStart:rowEnd)' * dz_mm;
    rowValues = log10( exp( -2*mu*zs ) ); 

    showImages = 1;
    if showImages == 1
      bscan_dB = intensity2dB( bscan );
      tmp_dB = shearImg( bscan_dB, shearAngle*pi/180 );
      figure; imshow( bscan_dB, [] );
      figure; imshow( tmp_dB, [min(bscan_dB(:)) max(bscan_dB(:))] );
      tmp_dB(1:rowStart-1,:)=0; tmp_dB(rowEnd+1:end,:)=0;
      tmp_dB(:,1:colStart-1)=0; tmp_dB(:,colEnd+1:end)=0;
      figure; imshow( tmp_dB, [min(bscan_dB(:)) max(bscan_dB(:))] );
    end

    figure; plotnice( rowValues, dataValues );
    Coefs = polyfit( rowValues, dataValues,1 );
    K = Coefs(1);
    disp(['K: ', num2str(K)]);
    hold on; plot( rowValues, K*rowValues+Coefs(2), 'k', 'LineWidth', 2 );

  end
end

