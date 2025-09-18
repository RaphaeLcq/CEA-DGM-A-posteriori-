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
*
*****************************************************************************/
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KsipB(i,j) = ∫_F ∇ψ_j · n_F * φ_i ds
% Mxy(i,j)  = ∫_F φ_i φ_j / |F|
% 
% SYNOPSIS [KsipB,Mxy]=KsipDGEdgeSimple(invh,oT,CoorBaryT,CoorEdgF,CoorMilF,EdgNormF)
%          
% INPUT  - invhT  : mise a l'echelle locale (inverse du diametre du triangle)
%        - oT : ordre d'approximation local au triangle
%        - CoorBaryT(1,2): coordonnees du barycentre du triangle
%        - CoorEdgF(2,2) : coordonnees de l'arete
%        - CoorMilF(1,2) : coordonnees du milieu de l'arete
%        - EdgNormF(1,2) : vecteur normal sortant a l'arete
% OUTPUT - KsipB : matrice ∇ψ_j · n_F * φ_i
%        - Mxy  : matrice de masse face
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KsipB,Mxy]=KsipDGEdgeNeumann(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF,EdgNormF)

% Matrice de masse arete du bord
Mxy=MassDGEdgeB(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF);

%-----------------------
% MATRICE ∇ψ_j · n_F * φ_i
%-----------------------

ndof=3*oT;
KsipB= zeros(ndof,ndof);

% vecteur normal scaled
Svec2  = invhT*EdgNormF;

% polynomes d'ordre 0
% Pour j=1 (constante): ∇ψ_1 = 0, donc KsipB(i,1) = 0
% Pour j=2,3 (lineaires): ∇ψ_j · n_F * φ_1 (constante)
KsipB(1,2) = Svec2(1);  % ∇ψ_2 · n_F * φ_1
KsipB(1,3) = Svec2(2);  % ∇ψ_3 · n_F * φ_1

% polynomes d'ordre 1
% KsipB(2,j) = ∫_F ∇ψ_j · n_F * φ_2 ds
KsipB(2,2) = Mxy(1,2)*Svec2(1);  % ∇ψ_2 · n_F * φ_2
KsipB(2,3) = Mxy(1,2)*Svec2(2);  % ∇ψ_3 · n_F * φ_2
KsipB(3,2) = Mxy(1,3)*Svec2(1);  % ∇ψ_2 · n_F * φ_3
KsipB(3,3) = Mxy(1,3)*Svec2(2);  % ∇ψ_3 · n_F * φ_3

if (oT>1)
  % interactions ordre 0 avec ordre 2
  KsipB(1,4) = KsipB(2,2);  % même gradient, φ different
  KsipB(1,5) = KsipB(2,3);
  KsipB(1,6) = KsipB(3,3);
  
  % polynomes d'ordre 2
  KsipB(2,4) = Mxy(2,2)*Svec2(1) + Mxy(2,3)*Svec2(2);
  KsipB(2,5) = 3*Mxy(2,2)*Svec2(1);  % coefficient pour terme quadratique
  KsipB(2,6) = Mxy(2,3)*Svec2(1) + Mxy(3,3)*Svec2(2);
  
  KsipB(3,4) = Mxy(2,3)*Svec2(1) + Mxy(3,3)*Svec2(2);
  KsipB(3,5) = KsipB(2,4);  % par structure des polynomes
  KsipB(3,6) = 3*Mxy(3,3)*Svec2(2);  % coefficient pour terme quadratique
  
  % polynomes d'ordre 3 et plus
  KsipB(4,4) = Mxy(2,4)*Svec2(1) + Mxy(2,6)*Svec2(2);
  KsipB(4,5) = Mxy(2,5)*Svec2(1) + 3*Mxy(2,4)*Svec2(2);
  KsipB(4,6) = 3*Mxy(2,6)*Svec2(1) + Mxy(3,6)*Svec2(2);
  
  KsipB(5,4) = KsipB(4,5);  % structure polynomiale
  KsipB(5,5) = 2*Mxy(2,5)*Svec2(1);  % terme de degré supérieur
  KsipB(5,6) = KsipB(4,4);  % structure polynomiale
  
  KsipB(6,4) = KsipB(4,6);  % structure polynomiale  
  KsipB(6,5) = KsipB(5,6);  % structure polynomiale
  KsipB(6,6) = 2*Mxy(3,6)*Svec2(2);  % terme de degré supérieur
end

end