% This software is offered under the GNU General Public License 3.0.  It 
% is offered without any warranty expressed or implied, including the 
% implied warranties of merchantability or fitness for a particular 
% purpose.


function [bscans] = load_cvol(cvol_file)

  complex_bscans = load_cvol_file(cvol_file);
  bscans = 20*log10(abs(complex_bscans));    
  bscans = squeeze( bscans(1,:,:) );
end

 function [cvol,hdr] = load_cvol_file(file)
% function = load_vol_file(file)
% Reads volume from .vol file
% .vol file is binary file with hdr
% hdr: hdrSize (uint16) = size of hdr (including hdr )
%         volSize (uint16) = volume size (z,x,y);
%         pos (uint16) = position of upper left corner (z,x,y);
% volume: vol (uint8) = size is specified by hdr.volSize

    warning('off','MATLAB:printf:BadEscapeSequenceInFormat');
    
% Open file
    fid = fopen(file,'r');

% Read hdr
    hdr.hdrSize = fread(fid,1,'uint16');
    tmphdrData = fread(fid,hdr.hdrSize-1,'uint16');
    hdr.volSize = double(tmphdrData(1:3));
    hdr.pos = tmphdrData(4:6);

% Read volume 
    cvol = double(zeros(hdr.volSize(1),hdr.volSize(2),hdr.volSize(3)));

    % Get the real part of the data
    for n = 1:hdr.volSize(3)
        data = fread(fid,hdr.volSize(1)*hdr.volSize(2),'double');
        cvol(:,:,n) = reshape(data,hdr.volSize(1),hdr.volSize(2));
    end

    % Get the imagianry part of the data
	for n = 1:hdr.volSize(3)
        data = fread(fid,hdr.volSize(1)*hdr.volSize(2),'double');
        cvol(:,:,n) = cvol(:,:,n)+1i*reshape(data,hdr.volSize(1),hdr.volSize(2));
    end

    fclose(fid);
    
 end

