% This software is offered under the GNU General Public License 3.0.  It 
% is offered without any warranty expressed or implied, including the 
% implied warranties of merchantability or fitness for a particular 
% purpose.

% Written by Nicholas Dwork (ndwork@stanford.edu) and 
% Gennifer Smith (gsmith9@stanford.edu).


function mu = muFitDRC( I, dz, hf, noisePower )

  Hres = 1 ./ hf(:) .* ( I(:).*I(:) ./ (I(:).*I(:) + noisePower ) );
  term = I .* Hres;

  nTerm = numel(term);
  integral = zeros(nTerm,1);
  for m = 1:(nTerm - 1)
    integral(m) = sum( term(m:end-1) );
  end

  mu = (.5/dz)*(term./integral);
  %mu = 0.5/dz * log( 1 + term./integral );

  mu( ~isfinite(mu) ) = 0;
end

