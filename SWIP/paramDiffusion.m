function D = paramDiffusion(alpha)
  
global Lshape SquareHarmonic Neumann NeumannTop LshapeNeumann SquareHole SquareSinus SquareSWIP

nonSWIP = Lshape + SquareHarmonic + Neumann + NeumannTop + LshapeNeumann + SquareHole + SquareSinus;
if (nonSWIP == 1)
  D = 1;
elseif (SquareSWIP == 1)
  D = SWIP_alphaToCoeffs(alpha);
  global DiffusionConstant = D;
endif
