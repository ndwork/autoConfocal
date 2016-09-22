
function surfaceLocs = findSurfaceLocs( in )

  sIn = size( in );
  nX = prod( sIn(2:end) );
  reshaped = reshape( in, [sIn(1) nX] );

  surfaceLocs = zeros( nX, 1 );
  for i=1:nX
    line = reshaped(:,i);
    surfaceLocs(i) = findSurfaceLoc( line );
  end

  if ~ismatrix(in)
    surfaceLocs = reshape( surfaceLocs, sIn(2:end) );
  end
end
