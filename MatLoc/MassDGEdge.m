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
% MassDGEdge.m:
%
% DG-FEM (SIPG monomes P1 ou P2)
%  
% Matrice de stabilite, fonctions de base monomiales pour une face donnee
%            (depend du triangle dans lequel la face est consideree)
%
% $\phi_i$ base locale dans $T$ (test fonction)
% $\phi_j$ base locale dans $T$ (trial function),
% $\vec{n}_F$ oriente vers l'exterieur de $T$
%
% Mxy(i,j)=  \sum_p w(p)\phi_i(XY_p)\,\phi_j(XY_p)
% 
% SYNOPSIS Mxy=MassDGEdge(invhT,oT,CoorBaryT,CoorMilF,wp3,Mat3T,wp5,Mat5T) 
%          
% INPUT : - invhT  : mise a l'echelle locale (inverse du diametre du triangle T)
%         - oT : ordre d'approximation local au triangle T  
%         - CoorBaryT(1,2) : coordonnees du barycentre de T
%         - CoorMilF(1,2) : coordonnees du milieu de l'arete
%         - wp3 (1,3)  : poids associes aux 3 points de la formule de Simpson
%         - Mat3T(3,5) : matrice des valeurs de x,y,xy,x^2,y^2 aux 3 points de la formule de Simpson
%         - wp5(1,5)   : poids associes aux 5 points de la formule de Boo
%         - Mat5T(5,5) : matrice des valeurs de x,y,xy,x^2,y^2 aux 5 points de la formule de Boo
% OUTPUT: - Mxy : matrice de masse aux faces
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mxy=MassDGEdge(invhT,oT,CoorBaryT,CoorMilF,wp3,Mat3T,wp5,Mat5T) 
  %
  Mxy(1,1  )=1;                          % 1*1
  Mxy(1,2:3)=invhT*(CoorMilF-CoorBaryT); % 1*x, 1*y
  % Polynomes d'ordre <=3
  %
  Mxy(2,2)=wp3*Mat3T(:,4); % x*x
  Mxy(2,3)=wp3*Mat3T(:,3); % x*y
  Mxy(3,3)=wp3*Mat3T(:,5); % y*y
  % SYMETRIE
  Mxy(2:3,1)=Mxy(1,2:3)';
  Mxy(3,2)  =Mxy(2,3);
  if (oT>=2)
    %
    Mxy(1,4)=Mxy(2,3); % 1*xy
    Mxy(1,5)=Mxy(2,2); % 1*xx
    Mxy(1,6)=Mxy(3,3); % 1*yy
    Mxy(2,4)=wp3*(Mat3T(:,1).*Mat3T(:,3)); % x*xy
    Mxy(2,5)=wp3*(Mat3T(:,1).*Mat3T(:,4)); % x*xx
    Mxy(2,6)=wp3*(Mat3T(:,1).*Mat3T(:,5)); % x*yy
    Mxy(3,4)=Mxy(2,6);                     % y*xy
    Mxy(3,5)=Mxy(2,4);                     % y*xx
    Mxy(3,6)=wp3*(Mat3T(:,2).*Mat3T(:,5)); % y*yy
    % SYMETRIE
    Mxy(4:6,1:3)=Mxy(1:3,4:6)';
    % Polynomes d'ordre<=6
    Mxy(4,4)=wp5*(Mat5T(:,3).*Mat5T(:,3)); % xy*xy
    Mxy(4,5)=wp5*(Mat5T(:,3).*Mat5T(:,4)); % xy*xx
    Mxy(4,6)=wp5*(Mat5T(:,3).*Mat5T(:,5)); % xy*yy
    Mxy(5,5)=wp5*(Mat5T(:,4).*Mat5T(:,4)); % xx*xx
    Mxy(5,6)=Mxy(4,4);                     % xx*yy
    Mxy(6,6)=wp5*(Mat5T(:,5).*Mat5T(:,5)); % yy*yy
    % SYMETRIE
    Mxy(5:6,4)=Mxy(4,5:6)';
    Mxy(6,5)=Mxy(5,6);
  end 
end
