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
% OUTPUT: - Kxy12 : matrice de consistance et symetrie
%         - Mxy1 : matrice de stabilite face T1
%         - Mxy2 : matrice de stabilite face T2
%         - Mxy12 : matrice de stabilite, couplage entre T1 et T2
%           couplage entre T1 et T2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mxy1,Mxy2,Mxy12]=MassDGEdge12Poly(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF,EdgNormF,LgF)
    CoorBaryT1=CoorBaryT(1,:);
    CoorBaryT2=CoorBaryT(2,:);
    invLgF=1/LgF;
    invhT1=invhT(1,1);
    invhT2=invhT(2,1); 
  %
     ndof1=3*oT(1); 
     ndof2=3*oT(2);
     Mxy1=zeros(ndof1,ndof1);
     Mxy2=zeros(ndof2,ndof2);
     Mxy12=zeros(ndof2,ndof1);
     %integral exacte dans T1
%     printf('x1=%7.2e\n',x1);
%     SiSj=CoorEdgF(2,:)-CoorEdgF(1,:);
%     LeSiSj=norm(SiSj);
%     Test=dot(bT1,SiSj)/LeSiSj;% fi fr2
%     if ( abs(abs(Test)) < 1.e-10)
%       printf('dot(bT1,SiSj)/|SiSj|=%7.2e\n',Test);
%     end
     %
%     if abs(x1)>1.e-10
%       CoorEdgF
%         XT1
%       [xp1,wp1,lambdaedg,npedg,rhoedg,thetaedg]=IntEdg_Boo5([CoorEdgF(1,:);(CoorEdgF(1,:)+CoorEdgF(2,:))/2]);
%       wp1=0.5*LgF*wp1;
%       [xp2,wp2,lambdaedg,npedg,rhoedg,thetaedg]=IntEdg_Boo5([(CoorEdgF(1,:)+CoorEdgF(2,:))/2;CoorEdgF(2,:)]);
%       wp2=0.5*LgF*wp2;
%       x2=I1_ex-(wp1*(exp(invnu*(xp1-XT1)*bT1).*(xp1(:,1)-CoorBaryT1(1,1)))+wp2*(exp(invnu*(xp2-XT1)*bT1).*(xp2(:,1)-CoorBaryT1(1,1))));
%       dif=abs(x2)-abs(x1);
%       %if dif>0
%        disp('fi chi mno sa7!!!!!!!!!!!!!!!!!!!!!!!!!!')
%      %end
%     end
     
     I1_x = integrate_exp_linear([0,0], [0,0], [1,0], CoorBaryT1,CoorEdgF(1,:), CoorEdgF(2,:));
     I1_y = integrate_exp_linear([0,0], [0,0], [0,1], CoorBaryT1,CoorEdgF(1,:), CoorEdgF(2,:));
     I1_xy=  integrate_b1_b2([0,1], CoorBaryT1, [1,0], CoorBaryT1, CoorEdgF(1,:), CoorEdgF(2,:));
     %I1_xy-wp*((xp(:,1)-CoorBaryT1(1,1)).*(xp(:,2)-CoorBaryT1(1,2))) 
     I1_yy=  integrate_b1_b2([0,1], CoorBaryT1, [0,1], CoorBaryT1, CoorEdgF(1,:), CoorEdgF(2,:));
     I1_xx=  integrate_b1_b2([1,0], CoorBaryT1, [1,0], CoorBaryT1, CoorEdgF(1,:), CoorEdgF(2,:));
     %integral exacte dans T2
     I2_x = integrate_exp_linear([0,0], [0,0], [1,0], CoorBaryT2,CoorEdgF(1,:), CoorEdgF(2,:));
     I2_y = integrate_exp_linear([0,0], [0,0], [0,1], CoorBaryT2,CoorEdgF(1,:), CoorEdgF(2,:));
     I2_xy=  integrate_b1_b2([0,1], CoorBaryT2, [1,0], CoorBaryT2, CoorEdgF(1,:), CoorEdgF(2,:));
     I2_yy=  integrate_b1_b2([0,1], CoorBaryT2, [0,1], CoorBaryT2, CoorEdgF(1,:), CoorEdgF(2,:));
     I2_xx=  integrate_b1_b2([1,0], CoorBaryT2, [1,0], CoorBaryT2, CoorEdgF(1,:), CoorEdgF(2,:));
     %integrale exacte dans T12
     Ixx=integrate_b1_b2([1,0], CoorBaryT1, [1,0], CoorBaryT2, CoorEdgF(1,:), CoorEdgF(2,:));
     Iyy=integrate_b1_b2([0,1], CoorBaryT1, [0,1], CoorBaryT2, CoorEdgF(1,:), CoorEdgF(2,:));
     Ix2y1=integrate_b1_b2([0,1], CoorBaryT1, [1,0], CoorBaryT2, CoorEdgF(1,:), CoorEdgF(2,:));
     Ix1y2=integrate_b1_b2([1,0], CoorBaryT1, [0,1], CoorBaryT2, CoorEdgF(1,:), CoorEdgF(2,:));
     %T1
     Mxy1(1,1)=1;
     Mxy1(2,1)=invLgF*I1_x*invhT1;
     Mxy1(1,2)=Mxy1(2,1);
     Mxy1(3,1)=invLgF*I1_y*invhT1;
     Mxy1(1,3)=Mxy1(3,1);
     Mxy1(2,2)=invLgF*invhT1^2*I1_xx;
     Mxy1(3,2)=invLgF*invhT1^2*I1_xy;
     Mxy1(2,3)=Mxy1(3,2);
     Mxy1(3,3)=invhT1^2*invLgF*I1_yy;
     %T2
     Mxy2(1,1)=1;
     Mxy2(2,1)=invLgF*I2_x*invhT2;
     Mxy2(1,2)=Mxy2(2,1);
     Mxy2(3,1)=invLgF*I2_y*invhT2;
     Mxy2(1,3)=Mxy2(3,1);
     Mxy2(2,2)=invLgF*invhT2^2*I2_xx;
     Mxy2(3,2)=invLgF*invhT2^2*I2_xy;
     Mxy2(2,3)=Mxy2(3,2);
     Mxy2(3,3)=invhT2^2*invLgF*I2_yy;
     %T12
     Mxy12(1,:)=-Mxy1(1,:);
     Mxy12(:,1)=-Mxy2(:,1);
     Mxy12(2,2)=-invhT1*invhT2*invLgF*Ixx;
     Mxy12(3,3)=-invhT1*invhT2*invLgF*Iyy;
     Mxy12(2,3)=-invLgF*invhT1*invhT2*Ix2y1;
     Mxy12(3,2)=-invLgF*invhT1*invhT2*Ix1y2;
end % end_function