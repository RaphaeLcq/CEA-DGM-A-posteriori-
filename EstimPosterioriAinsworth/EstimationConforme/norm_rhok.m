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
%% % Calcul des quantités faisant intervenir rho_K dans l'article p1794 eq (53)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Raphael Lecoq CEA
%
%
% SYNOPSIS [norm2_rhok, norm2_curl_rhok] = norm_rhok(idof,Uh)
%
% GLOBALS:
%   - ordre                      : ordre du schéma (P1 ou P2) 
%   - Nbtri                      : nombre de triangle
%   - TriEdg(Nbtri,3)            : indices globaux des faces du triangle
%   - EdgTri(Nbedg,2)            : indices globaux des triangles qui partagent une face
%   - Aires(Nbtri,1)             : aire d'un triangle
%   - EdgNorm(Nbtri,2)           : vecteur normal non normalisé, orienté vers l'extérieur du triangle
%   - LgEdg2(Nbedg,1)            : carré de la longueur de la face
% INPUT : 
%   - Uh                         : solution numérique 
%   - idof(Nbtri,3)              : index des degrés de liberté
% 
% OUTPUT:
%   - norm2_curl_rhok(Nbtri,1)   : C*_t^2 * Aire * carré de la norme L2 de curl(rho_T)
%   - norm2_rhok(Nbtri,1)        : carré de la norme L2 de rho_T 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [norm2_rhok, norm2_curl_rhok] = norm_rhok(idof,Uh)
  
global ordre 
global Nbtri TriEdg EdgTri 
global Aires EdgNorm LgEdg2 

rho_K = representerK(idof, Uh); %% [Nbtri x 3 x 2] coefficients p^K_n of eq (41) p 1792 without dividing by 2|K|

norm2_rhok = normRhok2(rho_K); %% [Nbtri,1]

[curlNorm2,Ck2] = computeCurlNorm2(rho_K); %% [Nbtri,1]

norm2_curl_rhok = Ck2.*curlNorm2.*Aires;


%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DEBUG
if (DEBUG)
  
  tolerance = 1e-5;
  [erreur_triangles, erreur_valeurs] = VisualiseDEBUGTriangles(idof, Uh, tolerance);

  Ind = find( (norm2_rhok -  norm2_curl_rhok) < 0);
  if isempty(Ind) == 0
    disp("%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%")
    disp("Indices où il y a un prb pour le calcul du rotationnel:\n")
    disp(Ind')
    disp("DIfférence, norme et curl norm pour le 1er indice")
    disp( norm2_rhok(Ind(1)) -  curlNorm2(Ind(1)) )
    disp(norm2_rhok(Ind(1)))
    disp(curlNorm2(Ind(1)))
    disp("%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%")
  endif
endif

endfunction
