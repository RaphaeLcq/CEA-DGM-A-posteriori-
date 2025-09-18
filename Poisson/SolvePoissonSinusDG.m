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
% SolvePoissonSinusDG.m:
%
% Convergence en maillage, EF DG P1 ou P2, base de monomes ou de Lagrange
%   Pb de Laplace dans un carre [0,1]*[0,1] :
%   -Delta Phi =(2*npi^2)*sin(npi*x)*sin(npi*y)
%   Solution Phi(x,y)=sin(npi*x)*sin(npi*y)
%
% SYNOPSIS [Eu0,Eu1,fig]=SolvePoissonSinusDG(fig)
%
% GLOBAL - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets (noeuds P1)
%        - CoorNeu2(Nbpt+Nbedg,2) : coordonnees (x, y) des noeuds P2
%        - RefNeu(Nbpt,1) : reference des sommets
%        - CoorBary(Nbtri,3) :coordonnees (x, y) des barycentres des triangles
%        - CoorMil(Nbedg,2)   : Coordonnees des milieux d'aretes
%		     - RefEdg(Nbedg,1) : Reference de chaque arete 
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets) 
%        - NumTri2(4*Nbtri,3) : liste de triangles du maillage P2
%                   (3 numeros de sommets)
%		     - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - EdgTri(Nbedg,2) : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere
%        - LgEdg2(Nbedg,1) : longueurs des aretes au carre
%        - EdgNorm(Nbedg,2) : vecteurs face-normale, orientes tri1->tri2
%        - Aires(Nbtri,1) : aires des triangles
%        - fig : numero de figure si visualisation
%        - ordre  : ordre d'approximation
%        - mi  : numero du maillage
% OUTPUT - Er0 : Erreur L2 normalisee \|\phi_h-\Pi_h(\phi)\|_0 /\|\Pi_h(\phi)\|_h
%        - Er1 : Erreur H1 normalisee   |\phi_h-\Pi_h(\phi))|_1/\|\Pi_h(\phi)\|_h
%        - fig = numero de la derniere figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Er0,Er1,fig]=SolvePoissonSinusDG(fig);
%
global CoorNeu CoorNeu2 RefNeu
global Nbtri CoorBary Aires NumTri NumTri2 TriEdg invDiaTri sigTri etaTri
global Nbedg CoorMil RefEdg LgEdg2 EdgNorm EdgTri NumEdg
global ordre mi 
global etaEdg
%
OrdTri=ordre*ones(Nbtri,1);
%
% degres de liberte pression et vitesse
% (dof=degrees of freedom)
ndof=3*OrdTri;
idof=zeros(Nbtri,3);
idof(:,3)=ndof; ndof_tot=sum(ndof);
idof(1,1)=1; idof(1,2)=ndof(1);
for t=2:Nbtri
  idof(t,1)=idof(t-1,2)+1;
  idof(t,2)=idof(t-1,2)+ndof(t);
end
%%%%%%%%%%%%%%%%%%%%%%%%
% renumerotation
%%%%%%%%%%%%%%%%%%%%%%%%
Profil=sparse(Nbtri,Nbtri);
for e=1:Nbedg
  tri=EdgTri(e,:);
  if (tri(2)>0)
    Profil(tri,tri)+=ones(2,2);
  end
end
%
renumT = symrcm(Profil);
sU=zeros(ndof_tot,1);
ideb=1;
for tt=1:Nbtri
  t=renumT(tt);
  i1=idof(t,1); ndoft=idof(t,3); 
  for p=0:ndoft-1
    sU(ideb+p,1)=i1+p;
  end
  ideb+=ndoft;
end
%%%%%%%%%%%%%%%%%%%%%%%%
% parametre de stabilite
%%%%%%%%%%%%%%%%%%%%%%%%
etaEdg=EtaParam(OrdTri);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrices pour la resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' Construction des matrices.\n');
[Ku,Mu]=KsipPolyDGsip(OrdTri,idof);
invMu=sparse(ndof_tot,ndof_tot);
for t=1:Nbtri
  debT=idof(t,1); finT=idof(t,2);
  % inversion locale (parallelisable)
  invMu(debT:finT,debT:finT)=inv(Mu(debT:finT,debT:finT));      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second membre  et solution exacte
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' Construction du second membre.\n');
[Phiexh,NormPhiex2h,NormGPhiex2h,RHS_Sinus]=SolutionPoissonSinusDG(Mu,invMu,Ku,idof);
VERIF=0;
if (VERIF)
  fprintf(' Verification du systeme lineaire.\n');
  errNum=abs((Ku*Phiexh-RHS_Sinus)'*Phiexh)/NormGPhiex2h;
  fprintf('P%i DG mesh_%i, errNum = %7.2e\n',ordre,mi,errNum);
  fprintf('------------------------------------------------\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resolution du systeme lineaire
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' Resolution.\n');
id = tic;
Ku_s=Ku(sU,sU); RHS_s=RHS_Sinus(sU,1);
%
% POUR VISUALISER L'INTERET DE LA RENUMEROTATION
if (fig>0)
  figure(fig);
  spy(Ku);
  title('Avant renumernotation');
  fig+=1;
  %
  figure(fig);
  spy(Ku_s);
  title('Apres renumernotation');
  fig+=1;
end
GCP=0; CH=0;
% Preconditionnement diagonal par blocs
% IDEE FUTURE : initialisation P^{k-1}_{dg} ou P^{k-1}_{lg}
if (GCP==1)
  invKu_s=sparse(ndof_tot,ndof_tot);
  for t=1:Nbtri
    i1=idof(t,1); i2=idof(t,2);
    i1_s=sU(i1); i2_s=sU(i2);
    Kuloc=Ku(i1:i2,i1:i2);
    invKu_s(i1_s:i2_s,i1_s:i2_s)=inv(Kuloc);
  end
  %
  [Phih_s,res,nit]=GCPDG(invKu_s,Ku_s,RHS_s);;
  fprintf('Nb it GCP = %i.\n',nit);
else
  if (CH==1)
    Kchol=chol(Ku_s);
    Phih_s=Kchol\(Kchol'\RHS_s);
  else
    Phih_s=Ku_s\RHS_s;
  end
end
%
Phih(sU,1)=Phih_s;
elapsed_time=toc(id);
fprintf('Temps de resolution (Ku*Phi=RHS) = %7.2e s\n',elapsed_time);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul des erreurs approche
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dPhi=Phih-Phiexh;
Er0=sqrt((dPhi'*Mu*dPhi)/NormGPhiex2h);
Er1=sqrt((dPhi'*Ku*dPhi)/NormGPhiex2h);
fprintf('P%i DG mesh_%i, ||Ph(Phiex)-Phih||_0/||Phiex||_1 = %7.2e\n',ordre,mi,Er0); 
fprintf('P%i DG mesh_%i, ||Ph(Phiex)-Phih||_h/||Phiex||_1 = %7.2e\n',ordre,mi,Er1);
fprintf('------------------------------------------------\n'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (fig>0)
  fig+=1;
  %% Matrice de masse LG
  nuT=ones(Nbtri,1);
  MLG = MatMassLG(ordre);
  Phih_LG=DGXYtoLG(Phih,ordre,OrdTri,MLG,idof);
  %
  [Phiex_LG] = SolutionPoissonSinusLG(ordre,MLG);
  %
  if (ordre==1)
    NT=NumTri; CN=CoorNeu;
  end
  if (ordre==2)
    NT=NumTri2; CN=CoorNeu2; 
  end
  tex=sprintf('Solution exacte Phi, P%i, mesh%i', ordre,mi);
  tDG=sprintf('Solution P%i DG Phi, mesh%i',ordre,mi);
  figure(fig)
  subplot(1,2,1)
  colormap ("jet");
  trisurf(NT,CN(:,1),CN(:,2),Phiex_LG);
  view(2);
  shading interp
  title(tex)
  colorbar;
  %
  subplot(1,2,2)
  trisurf(NT,CN(:,1),CN(:,2),Phih_LG);
  view(2);
  shading interp
  title(tDG)
  colorbar;
  %
end
