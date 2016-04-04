function image = loadLengBscan(fileName)

  height = 1024;
  width = 1024;
  nScans = 5;
  fid = fopen(fileName, 'rb');
  temp = fread(fid, inf, 'uint8');
  image = reshape(temp, [height,width,nScans]);

%  figure;
%  imshow(imrotate(squeeze(image(:,:,1)),90), []);
%   for n = 1:nScans
%     subplot(2,3,n); 
%     imshow(imrotate(squeeze(image(:,:,n)),90), []);
%     hold on;
%   end

end
