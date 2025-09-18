%{
/****************************************************************************
* Copyright (c) 2022, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Raphaël Lecoq CEA
%
%
% SYNOPSIS - Renvoie le nom du fichier de maillage pour le problème considéré 
%

% GLOBAL - Lshape SquareHarmonic Neumann NeumannTop LshapeNeumann SquareHole SquareSinus SquareSWIP  : problème considéré


function filename = NomFichierMaillage(mii)
  
global raffinement SquareSinus Lshape SquareHarmonic Neumann NeumannTop LshapeNeumann SquareHole SquareSWIP

if (raffinement == 0 && SquareSinus == 1)
  fprintf(' Maillage Square_h%i\n',mii);
  filename = sprintf('Square_h%i',mii);
elseif (raffinement == 1 && SquareSinus == 1)
  fprintf('Maillage SquareRaff_h%i\n',mii);
  filename = sprintf('SquareRaff%i',mii);
elseif (raffinement == 0 && Lshape == 1)
  fprintf(' Maillage Lshape non raffiné %i\n',mii);
  filename = sprintf('LshapeNonRaff%i',mii);
elseif (raffinement == 1 && Lshape == 1)
  fprintf(' Maillage Lshape%i\n',mii);
  filename = sprintf('Lshape%i',mii);
elseif (raffinement == 0 && SquareHarmonic == 1)
  fprintf(' Maillage Lshape%i\n',mii);
  filename = sprintf('Square_h%i',mii);
elseif (raffinement == 1 && SquareHarmonic == 1)
  fprintf(' Maillage Lshape%i\n',mii);
  filename = sprintf('SquareHarmonicTest%i',mii);
elseif (raffinement == 0 && Neumann == 1)
  fprintf('Maillage carré Neumann %i\n', mii);
  filename = sprintf('SquareNeumann_h%i',mii);
elseif (raffinement == 0 && NeumannTop == 1)
  fprintf('Maillage carré avec Neumann au top, dirichlet ailleurs %i\n', mii);
  filename = sprintf('SquareMixed_h%i',mii);
elseif (raffinement == 0 && LshapeNeumann == 1)
  fprintf('Maillage Lshape avec condition de Neumann %i\n', mii);
  filename = sprintf('LshapeNeumann_h%i',mii);
elseif (raffinement == 1 && LshapeNeumann == 1)
  fprintf('Maillage Lshape raffiné avec condition de Neumann %i\n', mii);
  filename = sprintf('LshapeNeumannRaff_h%i',mii);
elseif (raffinement == 0 && SquareHole == 1)
  fprintf('Maillage SquareHole %i\n', mii);
  filename = sprintf('SquareHole_h%i',mii);
elseif (raffinement == 1 && SquareHole == 1)
  fprintf('Maillage SquareHoleRaff %i\n', mii);
  filename = sprintf('SquareHoleRaff_h%i',mii);
elseif (raffinement == 0 && SquareSWIP == 1)
  fprintf('Maillage SquareSWIP %i\n',mii);
  filename = sprintf('SquareSWIP_h%i',mii);
elseif (raffinement == 1 && SquareSWIP == 1)
  fprintf('Maillage SquareSWIP Raffiné %i\n',mii);
  filename = sprintf('SquareSWIP_Raff_h%i',mii);
endif

endfunction