
function GradError2 = GradL2Error(Uh, idof)

global SquareSinus Neumann Lshape LshapeNeumann SquareHarmonic NeumannTop SquareHole SquareSWIP


disp('test')
if (SquareSinus == 1)
  GradError2 = GradL2ErrorSinus(Uh,idof);
elseif (Neumann == 1)
  GradError GradL2ErrorSquareNeumann(Uh, idof);
elseif (Lshape == 1)
  GradError2 = GradL2ErrorLshape(Uh,idof);
elseif (LshapeNeumann == 1)
  GradError2 = GradL2ErrorLshapeNeumann(Uh,idof);
elseif (SquareHarmonic == 1)
  GradError2 = GradL2ErrorHarmonicSquare(Uh, idof);
elseif (Neumann == 1)
  GradError2 = GradL2ErrorSquareNeumann(Uh,idof);
elseif (NeumannTop == 1)
  GradError2 = GradL2ErrorSquareTopNeumann(Uh,idof);
elseif (SquareHole == 1)
  disp('This problem is not yet implemented');
elseif (SquareSWIP == 1)
  GradError2 = GradL2ErrorSWIP_Square(Uh,idof);
endif

endfunction
