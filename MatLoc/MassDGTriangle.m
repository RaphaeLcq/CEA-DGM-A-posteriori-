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
% MassDGTriangle.m:
% Matrices de masse des monomes P1 et P2 pour un triangle donne
%
% SYNOPSIS Mass = MassDGTriangle(aire,invhT,oT,CoorNeuT,CoorMilT,CoorBaryT)
%
% Pour tester : aire=0.5; invh=1/sqrt(2); ordre=2;
%                         CoorNeu =[0,0;1,0;0,1]; 
%                         CoorMil =[0.5,0.5;0,0.5;0.5,0];
%                         CoorBary=[1/3,1/3];
%                         
%
% INPUT  - aire  : aire du triangle
%        - invhT  : mise a l'echelle (inverse du diametre du triangle)
%        - oT : ordre polynomial local
%        - CoorNeuT(3,2) : coordonnees des sommets du triangle
%        - CoorMilT(3,2) : coordonnees des milieux des aretes du triangle
%        - CoorBaryT(1,2): coordonnees du barycentre du triangle
%       
% OUTPUT - Mass: matrice de masse locale fonctions de base monomiales
%                [(x-xT)/hT]^nx*[(y-yT)/hT]^ny avec nx, ny<=ordre
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mass=MassDGTriangle(aire,invhT,oT,CoorNeuT,CoorMilT,CoorBaryT)
%
%% MATRICE DE MASSE
%
ndof=3*oT;
Mass=zeros(ndof,ndof);
%
% Polynome d'ordre 0
Mass(1,1)=aire;
% Polynomes d'ordre 1
  % on prend comme origine le barycentre du triangle
  % Mass(1,2:3)=invh*aire*(CoorBaryT-CoorOriT);
  %             NUL si CoorOriT=CoorBaryT
  % Mass(2:3,1)=Mass(1,2:3)';
% Polynomes d'ordre 2
% [xyp,wp,lambda,np]=IntTri_Ham3(XY);
wp=aire/3;
XMi=invhT*(CoorMilT(:,1)-ones(3,1)*CoorBaryT(1));
YMi=invhT*(CoorMilT(:,2)-ones(3,1)*CoorBaryT(2));
Mass(2,3)=wp*XMi'*YMi;
Mass(3,2)=Mass(2,3);
Mass(2,2)=wp*XMi'*XMi;
Mass(3,3)=wp*YMi'*YMi;
if (oT>1)
  % Ligne 1
  Mass(1,4)=Mass(2,3); % x*y
  Mass(1,5)=Mass(2,2); % x*x
  Mass(1,6)=Mass(3,3); % y*y
  % Polynomes d'ordre 3
  [xyp4,wp4,lambda4,np4]=IntTri_Ham4(CoorNeuT);
  xypT=invhT*(xyp4-ones(np4,1)*CoorBaryT);
  xypT2=xypT.^2;
  %  Ligne 2
  wp4=aire*wp4;
  Mass(2,4)=wp4*(xypT2(:,1).*xypT (:,2)); % x*x*y
  Mass(2,5)=wp4*(xypT2(:,1).*xypT (:,1)); % x*x*x
  Mass(2,6)=wp4*(xypT (:,1).*xypT2(:,2)); % x*y*y
  %  Ligne 3
  Mass(3,4)=Mass(2,6);      % x*y*y
  Mass(3,5)=Mass(2,4);      % x*x*y
  Mass(3,6)=wp4*(xypT2(:,2).*xypT(:,2)); % y*y*y
  % SYMETRIE
  Mass(4:6,1:3)=Mass(1:3,4:6)';
  %
  % Polynomes d'ordre 4
  [xyp7,wp7,lambda7,np7]=IntTri_Ham7(CoorNeuT);
  xypT=invhT*(xyp7-ones(np7,1)*CoorBaryT);
  xypT2=xypT.^2;
  xypT3=xypT.*xypT2;
  xypT4=xypT2.*xypT2;
  %
  wp7=aire*wp7;
  %  Ligne 4
  Mass(4,4)=wp7*(xypT2(:,1).*xypT2(:,2)); % x*x*y*y
  Mass(4,5)=wp7*(xypT3(:,1).*xypT (:,2)); % x*x*x*y
  Mass(4,6)=wp7*(xypT (:,1).*xypT3(:,2)); % x*y*y*y
  %
  Mass(5,5)=wp7*xypT4(:,1); % x*x*x*x
  Mass(5,6)=Mass(4,4);      % x*x*y*y
  %
  Mass(6,6)=wp7*xypT4(:,2); % y*y*y*y
  % SYMETRIE
  Mass(5:6,4)=Mass(4,5:6)';
  Mass(6  ,5)=Mass(5,6);
  
end