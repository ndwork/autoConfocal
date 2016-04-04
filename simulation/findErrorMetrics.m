
function meanEtbDepth = findErrorMetrics(muFit,trueMu,...
  z,offsetThreshPercent)

  [~, N] = size(muFit);

  etbDepths = zeros(1,N);
  vDepths = zeros(1,N);

  for i=1:N
    thisETB = findErrorEnergyTooBigDepth( offsetThreshPercent, ...
      trueMu(:,i), muFit(:,i), z );

    etbDepths(i) = thisETB;
  end

  meanEtbDepth = mean(etbDepths);
end

