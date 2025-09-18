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
% SYNOPSIS - Evalue fonction de Dirichlet du problème considéré en (x,y)
%

% GLOBAL - Lshape SquareHarmonic Neumann NeumannTop LshapeNeumann SquareHole SquareSinus SquareSWIP  : problème considéré
% Inputs - x(np,1) 
%        - y(np,1) 
% Output - fvals(np,1) : valeur de la donnée de Dirichlet en ce point


function value = eval_donneesDirichlet(x,y)
  
  global Lshape SquareHarmonic Neumann NeumannTop LshapeNeumann SquareHole SquareSinus SquareSWIP
  np = length(x);
  value = zeros(np,1);
  if (SquareSinus == 1)
    value = zeros(np,1); 
  elseif (Lshape == 1)
    [rp,thetap] = CartesianToPolarCentered(x,y);
    value = FonctionDiffusionLshape(rp,thetap);
  elseif (SquareHarmonic == 1)
    [rp,thetap] = CartesianToPolarCentered(x,y);
    value = FonctionHarmonicSquare(rp,thetap);
  elseif (Neumann == 1)
    value = zeros(np,1);
  elseif (NeumannTop == 1)
    value = zeros(np,1);
  elseif (LshapeNeumann == 1)
    value = zeros(np,1);
  elseif (SquareHole == 1);
    value = zeros(np,1);
  elseif (SquareSWIP == 1)
    [rp,thetap] = CartesianToPolarCentered(x,y);
    value = FonctionDiffusionSWIP(rp,thetap);
  endif

endfunction
