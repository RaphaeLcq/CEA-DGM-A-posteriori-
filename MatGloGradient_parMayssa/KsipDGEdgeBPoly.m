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
% KsipDGEdgeB.m:
% 
% DG-FEM (SIPG monomes P1 ou P2)
% 
% Matrice de consitance et stabilite
%
% $\phi_i$ base locale dans $T_1$ (test fonction)
% $\psi_j$ base locale dans $T_2$ (trial function)
%
% KsipB(i,j)=-0.5*\int_F\grad(\phi_i\,\phi_j)\cdot\nvec_F
% Mxy(i,j)  = \int_F \phi_i \phi_j / |F|
% 
% SYNOPSIS [KsipB,Mxy]=KsipDGEdgeB(invh,oT,CoorBaryT,CoorEdgF,CoorMilF,EdgNormF)
%          
% INPUT  - invhT  : mise a l'echelle locale (inverse du diametre du triangle)
%        - oT : ordre d'approximation local au triangle
%        - CoorBaryT(1,2): coordonnees du barycentre du triangle
%        - CoorEdgF(2,2) : coordonnees de l'arete
%        - CoorMilF(1,2) : coordonnees du milieu de l'arete
%        - EdgNormF(1,2) : vecteur normal sortant a l'arete
% OUTPUT - KxyB : matrice de consistance et symetrie
%        - Mxy  : matrice de stabilite face
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Kxy1,Mxy1]=KsipDGEdgeBPoly(invhT1,oT,CoorBaryT,CoorEdgF,CoorMilF,EdgNormF,LgF)
%-----------------------
% MATRICE DE CONSISTANCE
%-----------------------
%
  % Matrice de consistance et de symetrie
  Nvec=EdgNormF/LgF;
  n1=Nvec(1,1);
  n2=Nvec(1,2);
  invLgF=1/LgF;
  %
  ndof1=3*oT; 
  Mxy1=zeros(ndof1,ndof1);Kxy1=Mxy1;
 % Kxy_co=zeros(ndof1,ndof2); Kxy_sy=Kxy_co;
  %
  %Kxy_co=int_F {\grad \phi_i}[\phi_j] avec le jump définit de T1 à T2 dans ce code.
  %XT1=CoorBaryT;
  %YM=invnu*(xyp-ones(np,1)*CoorBaryT(1,2));
  %XT=[XM,YM];
  %T
  I1_x = integrate_exp_linear([0,0],[0,0], [1,0], CoorBaryT,CoorEdgF(1,:), CoorEdgF(2,:));
  I1_y = integrate_exp_linear([0,0],[0,0], [0,1], CoorBaryT,CoorEdgF(1,:), CoorEdgF(2,:));
  I1_xy=  integrate_b1_b2([0,1], CoorBaryT, [1,0], CoorBaryT, CoorEdgF(1,:), CoorEdgF(2,:));
  I1_yy=  integrate_b1_b2([0,1], CoorBaryT, [0,1], CoorBaryT, CoorEdgF(1,:), CoorEdgF(2,:));
  I1_xx=  integrate_b1_b2([1,0], CoorBaryT, [1,0], CoorBaryT, CoorEdgF(1,:), CoorEdgF(2,:));
  %K1(j,i)=0.5(int_F \grad phu_i.n phi_j+ \grad phu_j.n phi_i
  
  Kxy1(2,1)=n1*LgF*invhT1;
  Kxy1(3,1)=n2*LgF*invhT1;
  Kxy1(2,2)=2*n1*invhT1^2*I1_x;
  Kxy1(3,2)=invhT1^2*(n1*I1_y+n2*I1_x);
  Kxy1(3,3)=2*n2*invhT1^2*I1_y;
  Kxy1(1,2)=Kxy1(2,1);
  Kxy1(1,3)=Kxy1(3,1);
  Kxy1(2,3)=Kxy1(3,2);

  Mxy1(1,1)=1;
  Mxy1(2,1)=invLgF*I1_x*invhT1;
  Mxy1(1,2)=Mxy1(2,1);
  Mxy1(3,1)=invLgF*I1_y*invhT1;
  Mxy1(1,3)=Mxy1(3,1);
  Mxy1(2,2)=invLgF*invhT1^2*I1_xx;
  Mxy1(3,2)=invLgF*invhT1^2*I1_xy;
  Mxy1(2,3)=Mxy1(3,2);
  Mxy1(3,3)=invhT1^2*invLgF*I1_yy;
  end 