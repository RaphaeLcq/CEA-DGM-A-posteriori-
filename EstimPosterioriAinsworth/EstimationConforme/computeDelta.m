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
%% Calcul des Delta définis dans l'article p1792 eq (43)
% Coordonnées de rho_K \cdot N sur les sommets des faces du triangle
%
% SYNOPSIS delta = computeDelta(tri, idof, gk, Uh, AGLO, IGLO)
%
% GLOBALS:
%   - NumTri(Nbtri,3)    : liste des triangles (3 numéros de sommets)
%   - CoorNeu(Nbpt,2)    : coordonnées des sommets
%   - CoorBary(Nbtri,2)  : coordonnées des barycentres
%   - CoorMil(Nbedg,2)   : coordonnées des milieux des faces
%   - TriEdg(Nbtri,3)    : indices globaux des faces du triangle
%   - EdgTri(Nbedg,2)    : indices globaux des triangles qui partagent une face
%   - EdgUnit(Nbtri,2)   : vecteur normal normalisé, orienté vers l'extérieur du triangle
%   - invDiaTri(Nbtri,1) : inverse du diamètre du triangle
%   - mu(Nbtri,3)        : mu(tri,i) indique si la normale de face i du triangle est dirigée vers l'extérieur du triangle (== 1) ou l'intérieur (== -1)
%   - kappa(Nbtri,1) : diffusion dans le triangle
% INPUT : 
%   - tri                : indice du triangle 
%   - idof(Nbtri,3)      : index des degrés de liberté
%   - gk(3,3)            : valeur de la fonction de diffusion numérique dans le triangle
%   - Uh                 : solution numérique 
%   - EdgUnitT(3,2)      : vecteur normal normalisé, orienté vers l'extérieur du triangle, numéroté dans l'ordre local du triangle  
% OUTPUT:
%   - delta(3,3)         : coefficients rho_n^(K)  dans l'article 1792 eq (42)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function delta = computeDelta(tri, idof, gk, Uh, EdgUnitT)
  
global CoorNeu CoorBary CoorMil NumTri TriEdg EdgTri
global mu kappa
global EdgUnit invDiaTri

AGLO = TriEdg(tri,:);
IGLO = NumTri(tri,:);

%% Évaluation des delta aux sommets 
delta = zeros(3,3); %% delta(i,j) = delta_i(x_j)
debT = idof(tri,1);
finT = idof(tri,2);
Uh_T = Uh(debT:finT);
%



%% Loop to compute delta(i,j) = delta_i(x_j), i != j where i is the edge, j the vertice
for i_edge_loc = 1:3
  for j_vert_loc = 1:3
    if i_edge_loc!=j_vert_loc 
      j_vert_glo = IGLO(j_vert_loc);
      sigma = kappa(tri).*eval_grad_fct_base(CoorNeu(j_vert_glo,:),tri)' * Uh_T; %% value of sigma on opposite vertice j      size(EdgUnitT(i_edge_loc,:))
      sigmav_tri = dot(sigma,EdgUnitT(i_edge_loc,:));  %%   Def de sigma_{nu_K} p1780
      delta(i_edge_loc,j_vert_loc) =  gk(i_edge_loc,j_vert_loc) - sigmav_tri ; %% Article p1792 eq (43)
   endif
  endfor
endfor


endfunction