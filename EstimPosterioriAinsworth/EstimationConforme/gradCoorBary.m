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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Raphael Lecoq CEA
%
% 
% 
% SYNOPSIS gradVerticeUh = gradCoorBary(idof,U_vertices)
%
% Gradient d'une fonction P1 dans base LG 
%
% U_vertice restreint à T = sum_{i=1}^3 u_i \lambda_i dans base barycentrique P1
% grad(U_vertices) restreint à T = sum_{i=1}^3 u_i grad( \lambda_i )) = sum_{i=1}^3 u_i -(d |T|)^{-1} S_i avec S_i normale à la face opposée 
%
% GLOBALS:
%   - Nbtri            : nombre de triangles
%   - NumTri(Nbpt,3)   : indice global des sommets du triangle
%   - TriEdg(Nbpt,3)   : indice global des faces du triangle
%   - EdgTri(Nbtri,2)  : indice global des triangles partagés par la face
%   - Aires(Nbtri,1)   : aire du triangle
%   - EdgNorm(Nbedg,2) : vecteur normal non normalisé d'une face 
% INPUT : 
%   - idof(Nbtri,3)    : index des degrés de liberté
%   - U_vertices       : valeur d'une fonction P1 aux sommets
% 
% OUTPUT:
%   - gradVerticeUh    : estimateur d'erreur non conforme
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function gradVerticeUh = gradCoorBary(idof,U_vertices)

global Nbtri NumTri TriEdg EdgTri
global Aires EdgNorm 

gradVerticeUh = zeros(Nbtri,2);
d = 2; %% dimension 
for tri = 1:Nbtri
  %
    invdarea = 1/(d*Aires(tri)); %%  1/(d * |T|) dans la forme de grad( \lambda_i )
    IGLO = NumTri(tri,:); 
    AGLO = TriEdg(tri,:); 
    EdgNormT = EdgNorm(AGLO,:); 
    for iloc=1:3 %% Vecteur normal non normalisé mais orienté par le triangle
      if (EdgTri(AGLO(iloc),1)~=tri)
        EdgNormT(iloc,:)=-EdgNormT(iloc,:);
      end
    end    
    gradVerticeUh(tri,:) = -invdarea*sum( U_vertices(IGLO).*EdgNormT,1); % grad(Uh) = sum_{i=1}^3 u_i grad( \lambda_i )) = sum_{i=1}^3 -1/(d |T|)*u_i * S_i, Uh P1 
endfor

endfunction 