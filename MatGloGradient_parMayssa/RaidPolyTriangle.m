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

% Matrices de masse des monomes P1 et P2 pour un triangle donne
%
% SYNOPSIS Mass = RaidPolyTriangle(aire,invhT,oT,CoorNeuT,CoorMilT,CoorBaryT)
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
function K=RaidPolyTriangle(aire,invhT,oT,CoorNeuT,CoorMilT,CoorBaryT) %OK
  %
  %% MATRICE DE MASSE
  %
  X1=CoorNeuT(1,:);
  X2=CoorNeuT(2,:);
  X3=CoorNeuT(3,:);
  ndof=3*oT;
  K=zeros(ndof,ndof);
  if oT==0
    ndof=1;
    K=zeros(ndof,ndof);
    K(1,1)=0;
  else
    K=zeros(ndof,ndof);
    K(3,3)=invhT^2*aire;
    K(2,2)=invhT^2*aire;
  end
end























