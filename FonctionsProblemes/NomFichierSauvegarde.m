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
% SYNOPSIS - Renvoie le nom du fichier de sauvegarde des erreurs pour le problème considéré 
%

% GLOBAL - Lshape SquareHarmonic Neumann NeumannTop LshapeNeumann SquareHole SquareSinus SquareSWIP  : problème considéré


function savefilename = NomFichierSauvegarde(save_dir, m0, m1)
  
  global raffinement SquareSinus Lshape SquareHarmonic Neumann NeumannTop LshapeNeumann SquareHole SquareSWIP SWIP alpha lambda

if (SquareSinus== 1)
  savefilename = fullfile(save_dir, sprintf('donnees_erreurs_m%i_m%i_SquareLbd%i_Raff%i.csv',m0,m1,lambda,raffinement));
elseif (Lshape == 1)
  savefilename = fullfile(save_dir, sprintf('donnees_erreurs_m%i_m%i_Lshape_Raff%i.csv',m0,m1,raffinement));
elseif (SquareHarmonic == 1)
  global HarmonicSquare_alpha ;
  global HarmonicSquare_a;
  global HarmonicSquare_b;
  savefilename = fullfile(save_dir, sprintf('donnees_erreurs_m%i_m%i_SquareHarmonic_alpha%i_a=%i_b=%i_Raff%i.csv',m0,m1,HarmonicSquare_alpha,HarmonicSquare_a, HarmonicSquare_b, raffinement));
elseif (Neumann == 1)
  savefilename = fullfile(save_dir, sprintf('donnees_erreurs_m%i_m%i_Neumann_Raff%i.csv',m0,m1,raffinement));
elseif (NeumannTop == 1)
  savefilename = fullfile(save_dir, sprintf('donnees_erreurs_m%i_m%i_NeumannTop_Raff%i.csv',m0,m1,raffinement));
elseif (LshapeNeumann == 1)
  savefilename = fullfile(save_dir, sprintf('donnees_erreurs_m%i_m%i_LshapeNeumann_Raff%i.csv',m0,m1,raffinement));
elseif (SquareHole == 1)
  savefilename = fullfile(save_dir, sprintf('donnees_erreurs_m%i_m%i_SquareHole_Raff%i.csv',m0,m1,raffinement));
elseif (SquareSWIP == 1 && SWIP == 1)
  savefilename = fullfile(save_dir, sprintf('donnees_erreurs_m%i_m%i_SquareSWIP_alpha%i_Raff%i.csv',m0,m1,alpha,raffinement));
elseif (SquareSWIP == 1 && SWIP == 0)
  savefilename = fullfile(save_dir, sprintf('donnees_erreurs_m%i_m%i_SquareSIP_alpha%i_Raff%i.csv',m0,m1,alpha,raffinement));
endif

  
endfunction
