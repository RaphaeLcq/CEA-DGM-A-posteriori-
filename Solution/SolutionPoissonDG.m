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
% Author :  CEA
%
%
% SYNOPSIS [Phiexh,NormPhiex2h,NormGPhiex2h]=SolutionStokesSinusDG(Mu,invMu,Ku,idof)
% Construct the Right hand side of the Discontinuous Galerkin method.
%
% GLOBAL
%        - CoorNeu(Nbpt,2)    : coordinates (x, y) of vertices (nodes P1)
%        - Nbtri              : global number of triangles
%        - NumTri(Nbtri,3)    : list of triangles  (3 vertices number)
%        - CoorBary(Nbtri,2)  : coordinates of barycentres of triangles
%        - invDiaTri(Nbtri,2) : inverse of diametres of triangles
%        - Aires(Nbtri,1)     : area of triangles
%		     - TriEdg(Nbtri,3)    : for each triangle, TriEdg(l,i) is the number of vertice opposite from vertice NumTri(l,i)
%                              (3 number of edges - full matrix Nbtri x 3)
%		     - RefEdg(Nbedg,1)    : reference of each vertice
%        - NumEdg(NbEdg,2)    : number of 2 nodes of each vertice
%        - LgEdg(Nbedg,1)     : length of edges
%        - etaEdg(NbEdg,1)    : value of eta on each edge
%		     - EdgTri(Nbedg,2)    : for each vertice, EdgTri(a,:) gives the number of the 2 triangles sharing vertice
%                                 EdgTri(a,2) = 0 if a is on the domain border
%        - EdgNorm(Nbedg,2)   : vector normal-edge, oriented tri1->tri2
%        - EdgUnit(Nbedg,2)   : unitary vector normal-edge, oriented tri1->tri2
%        - invLgEdg(Nbedg,1)  : inverse of the edge length
%        - kappa(Nbtri,1)     : diffusion in the triangle
% INPUT  - Mu   : mass matrix
%        - invMu: inverse of mass matrix
%        - Ku   : stiffness matrix
%        - idof : indices of dofs (degrees of freedom)
% OUTPUT
%        - Phiexh :exact solution projected on polynomials
%        - NormPhiex2h : discrete norm L2
%        - NormGPhiex2h : broken gradient discrete norm
%        - RHS : Right hand side of the system
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Phiexh,NormPhiex2h,NormGPhiex2h,RHS]=SolutionPoissonDG(Mu,invMu,Ku,idof)
global CoorNeu  ordre
global Nbtri NumTri CoorBary invDiaTri Aires TriEdg RefEdg
global NumEdg LgEdg etaEdg  EdgTri EdgNorm EdgUnit invLgEdg
global kappa

%
ndof=idof(Nbtri,2);
RHS=zeros(ndof,1); Phiexh=zeros(ndof,1);

for t=1:Nbtri
  IGLO=NumTri(t,:); CoorNeuT=CoorNeu(IGLO,:);
  deb=idof(t,1); fin=idof(t,2); ndofT=idof(t,3);
  %
  % volume of the element
  aire=Aires(t);
  % Integration points
  [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
  awp=aire*wp';
  PhiF=evalDonnees_Source(xyp(:,1), xyp(:,2));
  %
  invhT=invDiaTri(t);
  CoorBaryT=ones(np,1)*CoorBary(t,:);
  xypT=invhT*(xyp-CoorBaryT);
  PsiDG=[ones(np,1),xypT];
  if (ndofT>3)
    PsiDG=[PsiDG,xypT(:,1).*xypT(:,2),xypT.^2];
  end
  RHS(deb:fin,1)+=sum(awp.*PhiF.*PsiDG,1)';
end

%%%%%%%%%%%%%% Bord dirichlet

%
EdgDiri=find(RefEdg>=1 & RefEdg <= 9)' ; %% Set of Dirichlet edges, Ref >= 10 are Neumann
EdgUnit=EdgNorm.*invLgEdg;

for edge=EdgDiri
  tri = EdgTri(edge,1); %%1er triangle car le 2eme triangle = 0 sur le bord
  debT = idof(tri,1);
  finT = idof(tri,2);

  IndVertices = NumEdg(edge,:);% Indices des deux extrémités de l’arête
  XY = CoorNeu(IndVertices,:);
  [xyp,wp,lambda,np]=IntEdg_Boo5(XY);
  awp = LgEdg(edge)*wp';
  #evaluation des fonctions de bases aux points de quadrature
  Phi_on_quad=eval_fct_base(xyp, tri); % taille : [ndof x nbr quadra]
  grad_phi = eval_grad_fct_base(xyp, tri); % taille : [ndof x 2 x nbr quadra]

  % evaluation des données du bord aux points de quadrature
  PhiDirichlet =  eval_donneesDirichlet(xyp(:,1), xyp(:,2));
  phi_times_data = Phi_on_quad' .* PhiDirichlet;
  phi_times_data =  phi_times_data;
  RHS(debT:finT) += etaEdg(edge)*etaKappa(edge) / LgEdg(edge) .* sum(awp .* phi_times_data)';

  % Gradients des fonctions de base au bord (évalués sur les pts de quadrature)
  grad_phi_dot_n = zeros(size(grad_phi,1), np);  % [ndof x nbr quadra]
  for i = 1:size(grad_phi,1)
    for k = 1:np
      grad_phi_dot_n(i,k) = dot(squeeze(grad_phi(i,:,k)), EdgUnit(edge,:)); %% squeeze pour transformer [1 x 2 x 1] ---> [2 x 1]
    endfor
  endfor
  % Terme d’intégrale \int g kappa * grad phi_i \cdot n
  grad_phi_times_data = etaKappa(edge) * grad_phi_dot_n' .* PhiDirichlet; % taille [ndof x np]
  RHS(debT:finT,1) -= sum( awp.*grad_phi_times_data)';
endfor


%%%%%%%%%%%%%% Bord Neumann

EdgNeum = find(RefEdg >= 10)'; %% Convention : Ref >=10 <=> Neumann

for edge = EdgNeum
  tri = EdgTri(edge,1);  %% Toujours le triangle intérieur
  debT = idof(tri,1);
  finT = idof(tri,2);
  % Points de quadrature sur l'arête
  IndVertices = NumEdg(edge,:);
  XY = CoorNeu(IndVertices,:);
  [xyp, wp, lambda, np] = IntEdg_Boo5(XY);
  awp = LgEdg(edge) * wp';

  PsiDG = eval_fct_base(xyp, tri);

  % Donnée Neumann g(x,y)
  g_Neumann = eval_donneesNeumann(xyp(:,1),xyp(:,2),edge);  % [np x 1]
  % Terme intégral Neumann : \int g * phi_i ds
  phi_times_g = PsiDG'.*g_Neumann;  % [np x ndof]
  weighted_phi_g =  awp.*phi_times_g ;  % [np x ndof]
  RHS(debT:finT,1) +=  sum(weighted_phi_g,1)';

endfor
##


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
  Phiexh(deb:fin,1)+= sum(PhiexT.*PsiDG,1)';
endfor

Phiexh = invMu*Phiexh;

%
NormGPhiex2h=Phiexh'*Ku*Phiexh;
NormPhiex2h=Phiexh'*Mu*Phiexh;

