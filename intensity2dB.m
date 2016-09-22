
function out = intensity2dB( in, direction )

  k = 20;

  if nargin < 2
    direction = 1;
  end

  if direction > 0
    out = k * log10( in );
  else
    out = 10.^( in / k );
  end

end
