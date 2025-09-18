
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
%% Calcul of the flux approximation gk Article p1790 eq (36)
% P1 approximation of the flux functions at the CoorMil
%
% SYNOPSIS gk = fluxFunctions(tri, idof, Uh)
%
% GLOBALS:
%   - NumTri(Nbtri,3)    : liste des triangles (3 numéros de sommets)
%   - CoorNeu(Nbpt,2)    : coordonnées des sommets
%   - CoorBary(Nbtri,2)  : coordonnées des barycentres
%   - etaEdg(Nbedg,1)    : valeur du paramètre de pénalisation eta sur la face
%   - TriEdg(Nbtri,3)    : indices globaux des faces du triangle
%   - EdgTri(Nbedg,2)    : indices globaux des triangles qui partagent une face
%   - EdgUnit(Nbtri,2)   : vecteur normal normalisé, orienté vers l'extérieur du triangle
%   - invLgEdg(Nbedg,1)  : inverse de la longueur de la face
%   - invDiaTri(Nbtri,1) : inverse du diamètre du triangle
%   - RefEdg(Nbedg,1)    : référence de Dirichlet/Neumann de la face 
%   - mu(Nbtri,3)        : mu(tri,i) indique si la normale de face i du triangle est dirigée vers l'extérieur du triangle (== 1) ou l'intérieur (== -1)
%   - kappa(Nbtri,1) : diffusion dans le triangle
% INPUT : 
%   - tri                : indice du triangle 
%   - Uh                 : solution numérique 
%   - idof(Nbtri,3)      : index des degrés de liberté
% 
% OUTPUT:
%   - rho_K(Nbtri,3,2)   : coefficients rho_n^(K)  dans l'article 1792 eq (42)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function gk = fluxFunctions(tri, idof, Uh)

global NumTri TriEdg EdgTri CoorNeu CoorBary  RefEdg
global etaEdg invDiaTri invLgEdg mu EdgUnit kappa 

global DEBUG

debT = idof(tri,1);
finT = idof(tri,2);
Uh_T = Uh(debT:finT);
IGLO = NumTri(tri,:);
AGLO = TriEdg(tri,:);
gk = zeros(3,3);

%
EdgUnitT=EdgUnit(AGLO,:);
for iloc=1:3
  if (EdgTri(AGLO(iloc),1)~=tri)
    EdgUnitT(iloc,:)=-EdgUnitT(iloc,:);
  endif
endfor
%

for i_edge_loc = 1:3
  for j_vertice_loc = 1:3
    if i_edge_loc != j_vertice_loc
      
      edg_glo = AGLO(i_edge_loc);
      vert_glo = IGLO(j_vertice_loc);
      
      if (RefEdg(edg_glo) < 10) %% Intérieur ou Dirichlet

        EdgeValue = eval_fct_base(CoorNeu(vert_glo,:),tri)' * Uh_T ;
        sigma = kappa(tri)*eval_grad_fct_base(CoorNeu(vert_glo,:),tri)' * Uh_T ;
        muk = mu(tri,i_edge_loc);
        [shared_tri, shared_EdgeValue, shared_sigma, shared_muk]  = shared_triangle(Uh, idof, tri, CoorNeu(vert_glo,:), edg_glo, vert_glo);
        jumpEdge = muk * EdgeValue + shared_muk * shared_EdgeValue;
        mean_sigma = sigma;
        if (RefEdg(edg_glo) == 0) %% no averaging on Dirichlet boundaries
          [avg1,avg2] = HarmonicAverageKappa(edg_glo);
          if (muk == 1)
            mean_sigma = avg1*sigma + avg2*shared_sigma;
          else
            mean_sigma = avg2*sigma + avg1*shared_sigma;
          endif
        endif
        oriented_EdgUnitT = muk .* EdgUnitT(i_edge_loc,:);
        mean_sigma_v = dot(mean_sigma,oriented_EdgUnitT );
        gk(i_edge_loc,j_vertice_loc) = muk * (mean_sigma_v - etaEdg(edg_glo) * etaKappa(edg_glo) * invLgEdg(edg_glo) * jumpEdge); %% eq (36) p 1790 /!\ muk = 1 sur les bords de Dirichlet 
        
      else %% Neumann
        gk(i_edge_loc,j_vertice_loc) = meanNeumannEdge(edg_glo);
      endif
      
    endif
  endfor
endfor

if (DEBUG)
  
endif 
endfunction