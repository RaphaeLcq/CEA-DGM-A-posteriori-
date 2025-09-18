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
% % Calcul de rho_K dans l'article 1792 eq (42)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Raphael Lecoq CEA
%
%
% SYNOPSIS rho_K = representerK(idof, Uh)
%
% GLOBALS:
%   - NumTri(Nbtri,3) : liste des triangles (3 numéros de sommets)
%   - Nbtri                      : nombre de triangle
%   - TriEdg(Nbtri,3)            : indices globaux des faces du triangle
%   - EdgTri(Nbedg,2)            : indices globaux des triangles qui partagent une face
%   - EdgUnit(Nbtri,2)           : vecteur normal normalisé, orienté vers l'extérieur du triangle
%   - LgEdg(Nbedg,1)             : longueur de la face
%   - RefEdg(Nbedg,1)            : référence de Dirichlet/Neumann de la face
%   - mu(Nbtri,3)                : mu(tri,i) indique si la normale de face i du triangle est dirigée vers l'extérieur du triangle (== 1) ou l'intérieur (== -1)
%   - kappa(Nbtri,1)         : diffusion dans le triangle

% INPUT :
%   - Uh                         : solution numérique
%   - idof(Nbtri,3)              : index des degrés de liberté
%
% OUTPUT:
%   - rho_K(Nbtri,3,2)           : coefficients rho_n^(K)  dans l'article 1792 eq (42)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function rho_K = representerK(idof, Uh)


global NumTri TriEdg EdgTri Nbtri RefEdg
global EdgNorm LgEdg EdgUnit
global mu kappa
%


mu = computeMu();
rho_K = zeros(Nbtri,3,2); %% 3 = nombre de faces, 2 = vecteur , rho_k(tri,:,:) coeffs de rho_k dans base barycentrique

for tri = 1:Nbtri
  IGLO = NumTri(tri,:); %% global index of vertices
  AGLO = TriEdg(tri,:); %% global index of opposite edges

  EdgNormT = EdgNorm(AGLO,:);
  for iloc=1:3
      if (EdgTri(AGLO(iloc),1)~=tri)
        EdgNormT(iloc,:)=-EdgNormT(iloc,:);
      endif
    endfor
  EdgUnitT=EdgUnit(AGLO,:);
  for iloc=1:3
      if (EdgTri(AGLO(iloc),1)~=tri)
        EdgUnitT(iloc,:)=-EdgUnitT(iloc,:);
      endif
    endfor
  EdgTanT = EdgNorm_to_EdgTan(tri, EdgNormT); %% Vecteurs tangentiels non normalisés
  gk = fluxFunctions(tri, idof, Uh);
  delta = computeDelta(tri, idof, gk, Uh, EdgUnitT);

  %% Coefficients barycentriques du représentant rho_K P1 p 1792 eq (42)
  %% delta(i,j) = delta_i(x_j)
  rho_K(tri,1,:) = LgEdg(AGLO(3)) .* delta(3,1) .* EdgTanT(2,:) - LgEdg(AGLO(2)) .* delta(2,1) .* EdgTanT(3,:);
  rho_K(tri,2,:) = LgEdg(AGLO(1)) .* delta(1,2) .* EdgTanT(3,:) - LgEdg(AGLO(3)) .* delta(3,2) .* EdgTanT(1,:);
  rho_K(tri,3,:) = LgEdg(AGLO(2)) .* delta(2,3) .* EdgTanT(1,:) - LgEdg(AGLO(1)) .* delta(1,3) .* EdgTanT(2,:);

endfor



%%%%%%%%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DEBUG
if(DEBUG) %% Vérification de gk + gk' = 0 preuve Lemme 5 p1790
  global Nbedg
  for edg=1:Nbedg
    if (RefEdg(edg)==0)
      tri_edg=EdgTri(edg,:);
      AGLO1 = TriEdg(tri_edg(1),:); %% global index of opposite edges in tri 1
      AGLO2 = TriEdg(tri_edg(2),:); %% global index of opposite edges in tri 2
      EdgUnitT1=EdgUnit(AGLO1,:);
      EdgUnitT2=EdgUnit(AGLO2,:);
      for iloc=1:3
        if (EdgTri(AGLO1(iloc),1)~=tri_edg(1))
          EdgUnitT1(iloc,:)=-EdgUnitT1(iloc,:);
        endif
      endfor
      for iloc=1:3
        if (EdgTri(AGLO2(iloc),1)~=tri_edg(2))
          EdgUnitT2(iloc,:)=-EdgUnitT2(iloc,:);
        endif
      endfor
      gk_t1 = fluxFunctions(tri_edg(1), idof, Uh, AGLO1, EdgUnitT1);
      gk_t2 = fluxFunctions(tri_edg(2), idof, Uh, AGLO2, EdgUnitT2);
      i_edg1=find(AGLO1==edg);
      i_edg2=find(AGLO2==edg);
        if(abs(sum(gk_t1(i_edg1,:))+sum(gk_t2(i_edg2,:)))>10e-6)
        disp("Error, flux are not equilibrated, see representerK")
      endif
    endif
  endfor


%%%%%%% Vérification de N_i \cdot rho_n^(K) = 2|K|Delta_i^(T)(x_n) dans preuve du Lemme 6 p1792
global Aires
for tri = 1:Nbtri

  DEBUG_IGLO = NumTri(tri,:); %% global index of vertices
  DEBUG_AGLO = TriEdg(tri,:); %% global index of opposite edges

  DEBUG_EdgNormT = EdgNorm(DEBUG_AGLO,:);
  for iloc=1:3
      if (EdgTri(DEBUG_AGLO(iloc),1)~=tri)
        DEBUG_EdgNormT(iloc,:)=-DEBUG_EdgNormT(iloc,:);
      endif
    endfor
  DEBUG_EdgUnitT=EdgUnit(DEBUG_AGLO,:);
  for iloc=1:3
      if (EdgTri(DEBUG_AGLO(iloc),1)~=tri)
        DEBUG_EdgUnitT(iloc,:)=-DEBUG_EdgUnitT(iloc,:);
      endif
    endfor
  EdgTanT = EdgNorm_to_EdgTan(tri, DEBUG_EdgNormT); %% Vecteurs tangentiels non normalisés
  gk = fluxFunctions(tri, idof, Uh, DEBUG_EdgUnitT);

  delta = computeDelta(tri, idof, gk, Uh, DEBUG_EdgUnitT);


  if ( abs( dot(DEBUG_EdgUnitT(1,:),squeeze(rho_K(tri,3,:))) - 2*Aires(tri)*delta(1,3) ) + abs( dot(DEBUG_EdgUnitT(1,:),squeeze(rho_K(tri,2,:))) - 2*Aires(tri)*delta(1,2) )>1e-10) ;
    disp("Probleme delta 2");
  endif
  if ( abs( dot(DEBUG_EdgUnitT(2,:),squeeze(rho_K(tri,3,:))) - 2*Aires(tri)*delta(2,3) ) + abs( dot(DEBUG_EdgUnitT(2,:),squeeze(rho_K(tri,1,:))) - 2*Aires(tri)*delta(2,1) )>1e-10 ) ;
    disp("Probleme delta 2");
  endif
  if ( abs( dot(DEBUG_EdgUnitT(3,:),squeeze(rho_K(tri,2,:))) - 2*Aires(tri)*delta(3,2) ) + abs( dot(DEBUG_EdgUnitT(3,:),squeeze(rho_K(tri,1,:))) - 2*Aires(tri)*delta(3,1) )>1e-10 ) ;
    disp("Probleme delta 3");
  endif

endfor

endif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

endfunction
