function fvals = evalSolution(x,y)

  global npi2 npi
  global Lshape SquareHarmonic Neumann NeumannTop LshapeNeumann SquareHole SquareSinus SquareSWIP
  np = length(x);
  fvals = zeros(np,1);
  if (SquareSinus == 1)
    fvals = sin(npi .* x) .* sin(npi .* y);
  elseif (Lshape == 1)
    [rp,thetap] = CartesianToPolarCentered(x,y);
    fvals = FonctionDiffusionLshape(rp,thetap);
  elseif (SquareHarmonic == 1)
    [rp,thetap] = CartesianToPolarCentered(x,y);
    fvals = FonctionHarmonicSquare(rp,thetap);
  elseif (Neumann == 1)
    fvals = FonctionSquareNeumann(x,y);
  elseif (NeumannTop == 1)
    fvals = FonctionSquareTopNeumann(x,y);
  elseif (LshapeNeumann == 1)
    [rp,thetap] = CartesianToPolarCentered(x,y);
    fvals = FonctionDiffusionLshapeNeumann(rp,thetap);
  elseif (SquareSWIP == 1);
    [rp,thetap] = CartesianToPolarCentered(x,y);
    fvals = FonctionDiffusionSWIP(rp,thetap);
  endif

endfunction
