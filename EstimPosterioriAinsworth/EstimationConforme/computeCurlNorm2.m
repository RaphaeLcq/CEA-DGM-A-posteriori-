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
% Author : Raphael Lecoq CEA
%
%% Computation of curl(rho_K) for eq (43) and (46) p 1793
%
% SYNOPSIS curlNorm2 = computeCurlNorm2(rho_K)
%
% GLOBALS:
%   - Nbtri                : nombre de triangle
%   - TriEdg(Nbtri,3)      : indices globaux des faces du triangle
%   - EdgTri(Nbedg,2)      : indices globaux des triangles qui partagent une face
%   - EdgNorm(Nbtri,2)     : vecteur normal non normalisé, orienté vers l'extérieur du triangle
%   - LgEdg(Nbedg,1)       : longueur de la face
%   - LgEdg2(Nbedg,1)      : carré de la longueur de la face
%   - ordre                : ordre du schéma (P1 ou P2)
%   - Aires(Nbtri,1)       : aire des triangles
%   - kappa(Nbtri,1)   : diffusion dans le triangle

% INPUT : 
%   - rho_K                : coefficients rho_n^(K)  dans l'article 1792 eq (42) 
%  
% OUTPUT:
%   - curlNorm2(Nbtri,3,2) : norme L2 au carré de Ck^2 * |T| *curl(rho_K) calculée avec l'article eq (43) and (46) p 1793
%   - Ck2(Nbtri,1)         : Constante C^*_t de l'article pour chaque triangle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [curlNorm2,Ck2] = computeCurlNorm2(rho_K)
  
global Nbtri TriEdg EdgTri 
global Aires LgEdg LgEdg2 EdgNorm
global ordre kappa
curlNorm2 = zeros(Nbtri,1);
Ck2 = zeros(Nbtri,1);
%
unDivisePar20 = 1/20; %% pour pas répéter division
for tri = 1:Nbtri

  area = Aires(tri);
  AGLO = TriEdg(tri,:); %% global index of opposite edges
  EdgNormT = EdgNorm(AGLO,:); 
  
  for iloc = 1:3
    if (EdgTri(AGLO(iloc),1) ~= tri)
      EdgNormT(iloc,:) = -EdgNormT(iloc,:);
    endif
  endfor
  curl_baryk = curlBary(area,EdgNormT);

  % Somme des produits scalaires rho_K^i ⋅ S_i / 2|T|
  dot_sum = sum(dot(squeeze(rho_K(tri, :, :))./ (2 * area),curl_baryk));
  % Norme L2 = aire * valeur^2
  curlNorm2(tri) = area * dot_sum^2;
  
  % Calcul de C^*_K(beta)^2 : aK/(20trace(Sk)) eq (20) p1783
  K = KriLGTriangle(area,EdgNormT,LgEdg2(TriEdg(tri,:)),ordre); %% p1781 eq (12)
  invTr = unDivisePar20/trace(K);
  ak = kappa(tri); %% Pour diffusion hétérogène
  Ck2(tri) = ak*invTr; 

  
endfor

endfunction
