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
%% Vecteurs tangentiels non normalisés dirigés dans le sens trigonométrique à partir des vecteurs normaux non normalisés
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : CEA
%
%
% SYNOPSIS EdgTan = EdgNorm_to_EdgTan(tri, EdgNormT) 
%
% GLOBALS:
%
% INPUT : 
%   - tri                         : indice du triangle 
%   - EdgNormT(3,2)              : vecteur normal non normalisé vers l'extérieur des faces du triangles
% 
% OUTPUT:
%   - EdgTan(3,2)  : vecteur tangent non normalisé des faces du triangles, orienté dans le sens trigonométrique par rapport aux normales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 

function EdgTan = EdgNorm_to_EdgTan(tri, EdgNormT) 

EdgTan = zeros(3,2);
EdgTan(1,:) = [-EdgNormT(1,2),EdgNormT(1,1)];
EdgTan(2,:) = [-EdgNormT(2,2),EdgNormT(2,1)];
EdgTan(3,:) = -EdgTan(2,:) - EdgTan(1,:);

   
endfunction