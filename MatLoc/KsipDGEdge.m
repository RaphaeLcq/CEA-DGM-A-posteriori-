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
%
% Author : Erell Jamelot CEA
%
% KsipDGEdge.m:
%
% DG-FEM (SIPG monomes P1 ou P2)
% 
% Matrices de stabilite T1 et 12
% Matrices de consistance, de symetrie et de stabilite
%           couplage entre T1 et T2
%
% $\phi_i$ base locale dans $T_1$ (test fonction)
% $\psi_j$ base locale dans $T_2$ (trial function)
%
% Kxy12(i,j)= -0.5*\int_F \grad\psi_j\cdot\nvec_F\,\phi_i  CONSISTANCE
%             +0.5*\int_F kappa1*omega1*\grad\phi_i\cdot\nvec_F\,\psi_j  SYMETRIE
% Mxy1(i,j) =  \int_F \phi_i \phi_j / |F| STABILITE
% Mxy2(i,j) =  \int_F \psi_i \psi_j / |F| STABILITE
% Mxy12(i,j)=  \int_F \phi_i \psi_j / |F|  STABILITE
% 
% SYNOPSIS [Kxy12,Mxy1,Mxy2,Mxy12]=KsipDGEdge(invhT,oT,CoorBaryT,...
%                                      ... CoorEdgF,CoorMilF,EdgNormF)
%          
% INPUT  - invhT(2)  : mise a l'echelle locale (inverse du diametre des triangles)
%        - oT(2) : ordre d'approximation local aux triangles T
%        - CoorBaryT(2,2): coordonnees des barycentres de triangles
%        - CoorEdgF(2,2) : coordonnees de l'arete
%        - CoorMilF(1,2) : coordonnees du milieu de l'arete
%        - EdgNormF(1,2)  : vecteur normal a l'arete
% OUTPUT - Kxy12 : matrice de consistance et symetrie
%        - Mxy1 : matrice de stabilite face T1
%        - Mxy2 : matrice de stabilite face T2
%        - Mxy12 : matrice de stabilite, couplage entre T1 et T2
%           couplage entre T1 et T2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Kxy12,Mxy1,Mxy2,Mxy12]=KsipDGEdge(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF,EdgNormF,kappa1,kappa2)
  %
  
  [Mxy1,Mxy2,Mxy12]=MassDGEdge12(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF);
  % Matrice de consistance et de symetrie
  

   kNvec  = 0.5*EdgNormF;
   Svec1  = kappa1*kNvec*invhT(1);
   Svec2  = -kappa2*kNvec*invhT(2);

  
  % 
  ndof1=3*oT(1); 
  ndof2=3*oT(2);
  Kxy_co=zeros(ndof1,ndof2); Kxy_sy=Kxy_co;
  %
  Kxy_co(:,2:3)= Mxy1(:,1)*Svec2;
  if (oT(2)>1)
    Kxy_co(:,4)  = Mxy12(:,3)*Svec2(1)+Mxy12(:,2)*Svec2(2);
    Kxy_co(:,5:6)= 2*Mxy12(:,2:3).*Svec2;
  end
  %
  Kxy_sy(2:3,:)= Svec1'*Mxy2(1,:);
  if (oT(1)>1)
    Kxy_sy(4,:)  = Mxy12(3,:)*Svec1(1)+Mxy12(2,:)*Svec1(2);
    Kxy_sy(5:6,:)= 2*Svec1'.*Mxy12(2:3,:);
  end
  %
  Kxy12=Kxy_co+Kxy_sy;

end
 

