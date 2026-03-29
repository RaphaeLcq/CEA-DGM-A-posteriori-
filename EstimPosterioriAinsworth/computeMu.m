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
% Author : CEA
%
%
% SYNOPSIS mu = computeMu()
%
% GLOBALS:
%   - TriEdg(Nbtri,3)    : indices globaux des faces du triangle
%   - EdgTri(Nbedg,2)    : indices globaux des triangles qui partagent une face
%   - Nbedg              : nombre de faces
%   - Nbtri              : nombre de triangles
% INPUT : 
%% OUTPUT:
%   - mu                 : mu(tri,i) indique si la normale de face i du triangle est dirigée vers l'extérieur du triangle (== 1) ou l'intérieur (== -1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mu = computeMu()
  
global RefEdg TriEdg EdgTri
global Nbedg Nbtri

mu = ones(Nbtri,3);

Ind = find(RefEdg == 0); %%On a mu = 1 sur le bord
for edge = Ind'
  
  tri = EdgTri(edge, 2); %% EdgTri(edge,2) est le + petit index global de EdgTri(edge,:) et les normales sont dirigées du triangle avec le + grand index vers le + petit dans EdgNorm
  loc = find(TriEdg(tri, :) == edge); %% "loc" indice local associé à la numérotation globale de la face "edge" dans le triangle "tri" considéré
  mu(tri,loc) = -1; 
endfor

  
endfunction
