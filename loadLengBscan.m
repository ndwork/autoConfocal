function image = loadLengBscan(fileName)

  height = 1024;
  width = 1024;
  nScans = 5;
  fid = fopen(fileName, 'rb');
  temp = fread(fid, inf, 'uint8');
  image = reshape(temp, [height,width,nScans]);

end
