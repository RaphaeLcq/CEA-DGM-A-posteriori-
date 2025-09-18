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
% MassDGEdge12.m:
%
% DG-FEM (SIPG monomes P1 ou P2)
% 
% Matrices de stabilite couplage, arete interieure
%
% $\phi_i$ base locale dans $T_1$ (test fonction)
% $\psi_j$ base locale dans $T_2$ (trial function)
%
% Mxy1(i,j) =  \int_F \phi_i \phi_j / |F| STABILITE
% Mxy2(i,j) =  \int_F \psi_i \psi_j / |F| STABILITE
% Mxy12(i,j)=  \int_F \phi_i \psi_j / |F|  STABILITE
% 
% SYNOPSIS [Mxy1,Mxy2,Mxy12]=MassDGEdge12(invhT,oT,CoorBaryT,...
%                                      ... CoorEdgF,CoorMilF,EdgNormF,)
%          
% INPUT : - invhT(2)  : mise a l'echelle locale (inverse du diametre des triangles)
%         - oT(2) : ordre d'approximation local aux triangles T
%         - CoorBaryT(2,2): coordonnees des barycentres de triangles
%         - CoorEdgF(2,2) : coordonnees de l'arete
%         - CoorMilF(1,2) : coordonnees du milieu de l'arete
% OUTPUT: - Mxy1 : matrice de stabilite face T1
%         - Mxy2 : matrice de stabilite face T2
%         - Mxy12 : matrice de stabilite, couplage entre T1 et T2
%           couplage entre T1 et T2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mxy1,Mxy2,Mxy12]=MassDGEdge12(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF)
 %
 [wp3,Mat3T1,Mat3T2,wp5,Mat5T1,Mat5T2]=ValPolyEdge(invhT,oT,CoorBaryT,CoorEdgF);
  %
  % On calcule les matrices de masses divisee par |F|
  %
  Mxy1=MassDGEdge(invhT(1),oT(1),CoorBaryT(1,:),CoorMilF,wp3,Mat3T1,wp5,Mat5T1) ;
  Mxy2=MassDGEdge(invhT(2),oT(2),CoorBaryT(2,:),CoorMilF,wp3,Mat3T2,wp5,Mat5T2) ;
  %
  ndof1=3*oT(1);ndof2=3*oT(2);
  %
  Mxy12=zeros(ndof1,ndof2);
  %
  Mxy12(1,1:3)=Mxy2(1,1:3); % 1*1,  1*x2, -*y2
  Mxy12(2:3,1)=Mxy1(2:3,1); % x1*1 ,y1*1
  % Polynomes d'ordre <=3
  Mxy12(2,2)=wp3*(Mat3T1(:,1).*Mat3T2(:,1)); % x1*x2
  Mxy12(2,3)=wp3*(Mat3T1(:,1).*Mat3T2(:,2)); % x1*y2
  Mxy12(3,2)=wp3*(Mat3T1(:,2).*Mat3T2(:,1)); % y1*x2
  Mxy12(3,3)=wp3*(Mat3T1(:,2).*Mat3T2(:,2)); % y1*y2
  o1o2=oT(1)+oT(2);
  if (o1o2>2)
    if (oT(1)>1)
      Mxy12(4:6,1)=Mxy1(4:6,1); % x1y1*1, x1x1*1, y1y1*1
      Mxy12(4,2)=wp3*(Mat3T1(:,3).*Mat3T2(:,1)); % x1y1*x2
      Mxy12(4,3)=wp3*(Mat3T1(:,3).*Mat3T2(:,2)); % x1y1*y2
      Mxy12(5,2)=wp3*(Mat3T1(:,4).*Mat3T2(:,1)); % x1x1*x2
      Mxy12(5,3)=wp3*(Mat3T1(:,4).*Mat3T2(:,2)); % x1x1*y2
      Mxy12(6,2)=wp3*(Mat3T1(:,5).*Mat3T2(:,1)); % y1y1*x2
      Mxy12(6,3)=wp3*(Mat3T1(:,5).*Mat3T2(:,2)); % y1y1*y2
    end
    if (oT(2)>1)
      Mxy12(1,4:6)=Mxy2(1,4:6); % 1*x2y2, 1*x2x2, 1*y2y2
      Mxy12(2,4)=wp3*(Mat3T1(:,1).*Mat3T2(:,3)); % x1*x2y2
      Mxy12(3,4)=wp3*(Mat3T1(:,2).*Mat3T2(:,3)); % y1*x2y2
      Mxy12(2,5)=wp3*(Mat3T1(:,1).*Mat3T2(:,4)); % x1*x2x2
      Mxy12(3,5)=wp3*(Mat3T1(:,2).*Mat3T2(:,4)); % y1*x2x2
      Mxy12(2,6)=wp3*(Mat3T1(:,1).*Mat3T2(:,5)); % x1*y2y2
      Mxy12(3,6)=wp3*(Mat3T1(:,2).*Mat3T2(:,5)); % y1*y2y2
    end
    % Polynomes d'ordre <=6
    if (o1o2>3)
      Mxy12(4,4)=wp5*(Mat5T1(:,3).*Mat5T2(:,3)); % x1y1*x2y2
      Mxy12(4,5)=wp5*(Mat5T1(:,3).*Mat5T2(:,4)); % x1y1*x2x2
      Mxy12(4,6)=wp5*(Mat5T1(:,3).*Mat5T2(:,5)); % x1y1*y2y2
      Mxy12(5,4)=wp5*(Mat5T1(:,4).*Mat5T2(:,3)); % x1x1*x2y2
      Mxy12(5,5)=wp5*(Mat5T1(:,4).*Mat5T2(:,4)); % x1x1*x2x2
      Mxy12(5,6)=wp5*(Mat5T1(:,4).*Mat5T2(:,5)); % x1x1*y2y2
      Mxy12(6,4)=wp5*(Mat5T1(:,5).*Mat5T2(:,3)); % y1y1*x2y2
      Mxy12(6,5)=wp5*(Mat5T1(:,5).*Mat5T2(:,4)); % y1y1*x2x2
      Mxy12(6,6)=wp5*(Mat5T1(:,5).*Mat5T2(:,5)); % y1y1*y2y2
    end
  end
end % end_function