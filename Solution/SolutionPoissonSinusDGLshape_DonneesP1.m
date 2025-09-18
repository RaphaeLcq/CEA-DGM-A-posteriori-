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
% Author : Raphael ecoq CEA
%
% RHS de DGM avec projection des données de Dirichlet dans les fonctions constantes par morceaux sur le bord
%
%
% SYNOPSIS [Phiexh,NormPhiex2h,NormGPhiex2h]=SolutionStokesSinusDG(Mu,invMu,Ku,idof)
%
% GLOBAL - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets (noeuds P1)
%        - CoorBary(Nbtri,2)   : Coordonnees des barycentres de triangles
%        - invDiaTri(Nbtri,2)   : inverse des diametres des triangles
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets) 
%		     - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%        - Aires(Nbtri,1) : aires des triangles
%        - ordre  : ordre d'approximation
%
% INPUT  - Mu   : matrice de masse
%        - invMu: inverse de la matrice de masse
%        - Ku   : matrice de raideur
%        - idof : indices des dof (degrees of freedom)
% OUTPUT - [Phiexh,NormPhiex2h,NormGPhiex2h,] : la solution exacte projetee sur les polynomes, les normes discretes et le second membre
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Phiexh,NormPhiex2h,NormGPhiex2h,RHS]=SolutionPoissonSinusDGLshape_DonneesP1(Mu,invMu,Ku,idof)
  
global DEBUG;

global Nbtri NumTri CoorBary invDiaTri Aires TriEdg RefEdg NumEdg LgEdg etaEdg  EdgTri EdgNorm EdgUnit invLgEdg
global DonneeDiri_P1 CoorNeu

%
EdgDiri=find(RefEdg>=1 & RefEdg <= 9)'; %% Set of Dirichlet edges
ndof=idof(Nbtri,2);
RHS=zeros(ndof,1); Phiexh=zeros(ndof,1);
EdgUnit=EdgNorm.*invLgEdg;


for t=1:Nbtri
  IGLO=NumTri(t,:); CoorNeuT=CoorNeu(IGLO,:);
  deb=idof(t,1); fin=idof(t,2); ndofT=idof(t,3);
  %   
  % volume of the element
  aire=Aires(t);
  % Integration points
  [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
  awp=aire*wp';
  PhiexT=awp.*evalDonnees_Source(xyp(:,1), xyp(:,2));
  %
  invhT=invDiaTri(t);
  CoorBaryT=ones(np,1)*CoorBary(t,:);
  xypT=invhT*(xyp-CoorBaryT);
  PsiDG=[ones(np,1),xypT];
  if (ndofT>3)
    PsiDG=[PsiDG,xypT(:,1).*xypT(:,2),xypT.^2];
  end
  RHS(deb:fin,1)+=sum(PhiexT.*PsiDG,1)';
end


for edge=EdgDiri
  
  IndVertices = NumEdg(edge,:);% Indices des deux extrémités de l'arête
  XY = CoorNeu(IndVertices,:); % 2x2 coordonnées des sommets
  %Evaluation des points de quadrature sur F
  [xyp,wp,lambda,np]=IntEdg_Sim3(XY);
  awp = LgEdg(edge)*wp';
  tri = EdgTri(edge,1); %%1er triangle car le 2eme triangle = 0 sur le bord
  debT = idof(tri,1);
  finT = idof(tri,2);
  
  % Evaluation des fonctions de bases aux points de quadrature
  Phi_on_quad = eval_fct_base(xyp, tri, CoorBary(tri,:), invDiaTri(tri)); % taille : [ndof x nbr quadra]
  grad_phi = eval_grad_fct_base(xyp, tri, CoorBary(tri,:), invDiaTri(tri)); % taille : [ndof x 2 x nbr quadra]
  
  % Calcul des valeurs de Dirichlet aux points de quadrature
  % DonneeDiri_P1 contient les coefficients P1 aux sommets globaux
  % Il faut récupérer les coefficients aux deux sommets de l'arête
  i = NumEdg(edge,1);  % Premier sommet global de l'arête
  j = NumEdg(edge,2);  % Deuxième sommet global de l'arête
  
  % Evaluation de g_D aux points de quadrature
  PhiDirichlet = zeros(np, 1);
  for k = 1:np
    PhiDirichlet(k) = DonneeDiri_P1(i) * lambda(k,1) + DonneeDiri_P1(j) * lambda(k,2);
  end
  
  % Premier terme : \int phi_i * g
  phi_times_data = Phi_on_quad' .* PhiDirichlet;  % taille [np x ndof]
  weighted_phi_times_data = awp .* phi_times_data; % awp est le vecteur des poids
  RHS(debT:finT) += etaKappa(edge) * etaEdg(edge) ./ LgEdg(edge) .* sum(weighted_phi_times_data)'; 
  
  % Gradients des fonctions de base au bord (évalués sur les pts de quadrature)
  grad_phi_dot_n = zeros(size(grad_phi,1), np);  % [ndof x nbr quadra]
  for i = 1:size(grad_phi,1)
    for k = 1:np
      grad_phi_dot_n(i,k) = dot(squeeze(grad_phi(i,:,k)), EdgUnit(edge,:)); %% squeeze pour transformer [1 x 2 x 1] ---> [2 x 1]
    endfor
  endfor
  
  % Terme d'intégrale \int g grad phi_i \cdot n
  grad_phi_times_data = etaKappa(edge) * grad_phi_dot_n' .* PhiDirichlet; % taille [np x ndof]
  RHS(debT:finT) -= sum(awp .* grad_phi_times_data)';
endfor


%%%%%%%%%%%%%% Bord Neumann

EdgNeum = find(RefEdg >= 10)'; %% Convention : Ref >=10 <=> Neumann

for edge = EdgNeum

  tri = EdgTri(edge,1);  %% Toujours le triangle à gauche
  debT = idof(tri,1);
  finT = idof(tri,2);

  % Points de quadrature sur l'arête
  IndVertices = NumEdg(edge,:);
  XY = CoorNeu(IndVertices,:);
  [xyp, wp, lambda, np] = IntEdg_Boo5(XY);
  awp = LgEdg(edge) * wp';

  % Fonction de base sur l'arête
  Phi_on_quad = eval_fct_base(xyp, tri, CoorBary(tri,:), invDiaTri(tri));  % [ndof x np]

  % Donnée Neumann g(x,y)
  g_Neumann = eval_donneesNeumann(xyp(:,1), xyp(:,2));  % [np x 1]

  % Terme intégral Neumann : \int g * phi_i ds
  phi_times_g = +Phi_on_quad'; % [np x ndof]
  weighted_phi_g = awp .* phi_times_g;     % [np x ndof]
  RHS(debT:finT) += sum(weighted_phi_g)';

endfor



%%% Phiexh = Projection sur DG de la solution exacte
for tri = 1:Nbtri
  IGLO=NumTri(tri,:); CoorNeuT=CoorNeu(IGLO,:);
  deb=idof(tri,1); fin=idof(tri,2); ndofT=idof(tri,3);
  %   
  % volume of the element
  aire=Aires(tri);
  % Integration points
  [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
  awp=aire*wp';
  PhiexT=awp.*evalSolution(xyp(:,1),xyp(:,2));
  %
  invhT=invDiaTri(tri);
  CoorBaryT=ones(np,1)*CoorBary(tri,:);
  xypT=invhT*(xyp-CoorBaryT);
  PsiDG=[ones(np,1),xypT];
  if (ndofT>3)
    PsiDG=[PsiDG,xypT(:,1).*xypT(:,2),xypT.^2];
  end
  Phiexh(deb:fin,1)+=sum(PhiexT.*PsiDG,1)';
endfor




Phiexh = invMu*Phiexh;
%
NormGPhiex2h=Phiexh'*Ku*Phiexh;
NormPhiex2h=Phiexh'*Mu*Phiexh;

