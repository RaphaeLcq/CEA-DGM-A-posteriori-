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
% Author : Mayssa Mroueh CEA
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
% Kxy12(i,j)= -0.5*\int_F\grad\psi_j\cdot\nvec_F\,\phi_i  CONSISTANCE
%             +0.5*\int_F\grad\phi_i\cdot\nvec_F\,\psi_j  SYMETRIE
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
function [Kxy1,Kxy2,Kxy12,Mxy1,Mxy2,Mxy12]=KsipDGEdgePoly(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF,EdgNormF,LgF,kappaAvg1,kappaAvg2,kappa1,kappa2)
  %
%  global NumTri CoorNeu
%  IGLO1=NumTri(T12(1,1),:);
%  CoorNeuT1=CoorNeu(IGLO1,:);
%  IGLO2=NumTri(T12(1,2),:);
%  CoorNeuT2=CoorNeu(IGLO2,:);
%  XS1=max(CoorNeuT1(:,1));
%  YS2=max(CoorNeuT1(:,2));
%  bST1=0.5*[invhT(1,1),invhT(1,1)];
  CoorBaryT1=CoorBaryT(1,:);
  CoorBaryT2=CoorBaryT(2,:);
  invhT1=invhT(1,1);
  invhT2=invhT(2,1);
  if oT==0
    ndof1=1; 
    ndof2=1;
    Kxy12=-1;
  else
  ndof1=3*oT(1); 
  ndof2=3*oT(2);
 end  
 Kxy1=zeros(ndof1,ndof1);
 Kxy2=zeros(ndof2,ndof2);
 Kxy12=zeros(ndof2,ndof1);
 % Matrice de consistance et de symetrie
 Nvec  = EdgNormF/LgF;
 n1=kappa1*kappaAvg1*Nvec(1,1);
 n2=kappa2*kappaAvg2*Nvec(1,2);
 %norm(2*Nvec)%=1 Ok
  %T1
  I1_x = integrate_exp_linear([0,0],[0,0], [1,0], CoorBaryT1,CoorEdgF(1,:), CoorEdgF(2,:));
  I1_y = integrate_exp_linear([0,0],[0,0], [0,1], CoorBaryT1,CoorEdgF(1,:), CoorEdgF(2,:));
  %K1(j,i)=0.5(int_F \grad phu_i.n phi_j+ \grad phu_j.n phi_i
  Kxy1(2,1)=n1*LgF*invhT1;
  Kxy1(3,1)=n2*LgF*invhT1;
  Kxy1(2,2)=2*n1*invhT1^2*I1_x;
  Kxy1(3,2)=invhT1^2*(n1*I1_y+n2*I1_x);
  Kxy1(3,3)=2*n2*invhT1^2*I1_y;
  Kxy1(1,2)=Kxy1(2,1);
  Kxy1(1,3)=Kxy1(3,1);
  Kxy1(2,3)=Kxy1(3,2);
  %T2
  %Kxy2(i,j)=-0.5*(\int_F \grad phi_j.nvec* phi_i+\grad phi_i.nvec* phi_j); u=phi_j; v=phi_i
  I2_x = integrate_exp_linear([0,0],[0,0], [1,0], CoorBaryT2,CoorEdgF(1,:), CoorEdgF(2,:));
  I2_y = integrate_exp_linear([0,0],[0,0], [0,1], CoorBaryT2,CoorEdgF(1,:), CoorEdgF(2,:));
  %K1(j,i)=0.5(int_F \grad phu_i.n phi_j+ \grad phu_j.n phi_i
  Kxy2(2,1)=-n1*LgF*invhT2;
  Kxy2(3,1)=-n2*LgF*invhT2;
  Kxy2(2,2)=-2*n1*invhT2^2*I2_x;
  Kxy2(3,2)=-invhT2^2*(n1*I2_y+n2*I2_x);
  Kxy2(3,3)=-2*n2*invhT2^2*I2_y;
  Kxy2(1,2)=Kxy2(2,1);
  Kxy2(1,3)=Kxy2(3,1);
  Kxy2(2,3)=Kxy2(3,2);
  %Couplage T1T2
  %Kxy1(i,j)=0.5*(-\int_F \grad phi_j.nvec* psi_i+\grad psi_i.nvec* phi_j); u=phi_j; v=phi_i
  Kxy12(1,:)=-Kxy1(1,:);%-wp*(dot(bT1,Nvec)*exp(invnu*(XT1*bT1))*invhT1*invnormbT1);
  Kxy12(:,1)=-Kxy2(:,1);
  Kxy12(2,2)=invhT1*invhT2*(-n1*I2_x+n1*I1_x);
  Kxy12(3,2)=invhT1*invhT2*(-n1*I2_y+n2*I1_x);
  Kxy12(2,3)=invhT1*invhT2*(-n2*I2_x+n1*I1_y);
  Kxy12(3,3)=invhT1*invhT2*(-n2*I2_y+n2*I1_y);
  %Matrice de Stabilité pour la vitesse:
  [Mxy1,Mxy2,Mxy12]=MassDGEdge12Poly(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF,EdgNormF,LgF);

end
 

