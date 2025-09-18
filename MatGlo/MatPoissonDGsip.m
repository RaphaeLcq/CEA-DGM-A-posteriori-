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
% MatPoissonDGsip.m:
%
% Matrices de raideur et matrice de masse
% Formulation DG SIP, fonctions de base monomiales
%
% SYNOPSIS [Ksip,Mass]=MatPoissonDGsip(OrdTri,idof)
%
% \phi_i et \phi_j de support T
%
% Ksip(i,j)= -0.5*\int_T(\Delta\phi_i\,\phi_j+\Delta\phi_j\,\phi_i)
%            -0.5*\sum_{F\in\pa T\cap\pa\Omega}\int_F\grad(\phi_i\,\phi_j)\cdot\nvec_F
%            +\sum_{F\in\pa T} eta_F\int_F\phi_i\,\phi_j
%
% \phi_i et \phi_j de support Ti et Tj, F=T_i\cap T_j, \nvec_F vecteur sortant de Ti
%
% Ksip(i,j)= 0.5*\int_F(\phi_j\grad\phi_i+\phi_i\grad\phi_j)\cdot\nvec_F
%            -eta_F\int_F\phi_i\,\phi_j

% OUTPUT - Ksip : matrice  de raideur, associee a l'operateur -Delta u, fonctions de base monomiales
%        - Mass : matrice de masse, fonctions de base monomiales
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ksip,Mass]=MatPoissonDGsip(OrdTri,idof)
%
global CoorNeu
global Nbtri CoorBary Aires NumTri TriEdg invDiaTri
global Nbedg NumEdg CoorMil RefEdg EdgNorm EdgTri LgEdg
global etaEdg 
%
Ndof=idof(Nbtri,2);
Ksip=sparse(Ndof,Ndof);
Mass=sparse(Ndof,Ndof);
%
%---------------------------------
% MATRICES D'INTEGRATION AUX FACES
%---------------------------------
%
for f=1:Nbedg
  ij=NumEdg(f,:);
  CoorEdgF=CoorNeu(ij,:);
  CoorMilF=CoorMil(f,:);
  EdgNormF=EdgNorm(f,:);
  T12=EdgTri(f,:);
  T1=T12(1); deb1=idof(T1,1); fin1=idof(T1,2); %ndof1=idof(T1,3);
  etaF=etaEdg(f); 
  etaK = etaKappa(f);
  %%[kappaAvg1, kappaAvg2] = HarmonicAverageKappa(f); On n'arrive pas encore à l'implémenter dans ce cas
  
  kappa1 = kappa(T1); 
  if (RefEdg(f)==0)
    % Face interieure
    % Kxy12(i,j)= -0.5 * kappa_2 * \int_F  \grad\psi_j\cdot\nvec_F\,\phi_i  CONSISTANCE
    %             +0.5 * kappa_1 * \int_F \grad\phi_i\cdot\nvec_F\,\psi_j  SYMETRIE
    % Mxy1(i,j) =  \int_F \phi_i \phi_j
    % Mxy2(i,j) =  \int_F \psi_i \psi_j
    % Mxy12(i,j)=  \int_F \phi_i \psi_j
    T2=T12(2); deb2=idof(T2,1); fin2=idof(T2,2);
    invhT=invDiaTri(T12);
    oT=OrdTri(T12);
    kappa2 = kappa(T2);
    CoorBaryT=CoorBary(T12,:);
    [Kxy12,Mxy1,Mxy2,Mxy12]=KsipDGEdge(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF,EdgNormF,kappa1,kappa2);
    %
    Ksip(deb1:fin1,deb1:fin1)+=(etaF*etaK*Mxy1);
    Ksip(deb2:fin2,deb2:fin2)+=(etaF*etaK*Mxy2);
    Ksip(deb1:fin1,deb2:fin2)+=(Kxy12-etaF*etaK*Mxy12);
    % SYMETRIE
    Ksip(deb2:fin2,deb1:fin1)=Ksip(deb1:fin1,deb2:fin2)';
  elseif (RefEdg(f) > 0 & RefEdg(f) < 10)
    % Face du bord Dirichlet
    % KsipB(i,j)=-0.5*\int_F kappa * \grad(\phi_i\,\phi_j)\cdot\nvec_F
    %            +etaF*kappa*\int_F\phi_i\,\phi_j/|F|
    invhT=invDiaTri(T1);
    oT=OrdTri(T1);
    CoorBaryT=CoorBary(T1,:);
    % Matrices de raideur et de masse arete du bord (consistance, symetrie et stabilite)
    [KsipB,Mxy]=KsipDGEdgeB(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF,EdgNormF);
    Ksip(deb1:fin1,deb1:fin1)+=(kappa1*KsipB+etaF*kappa1*Mxy);
  elseif (RefEdg(f) >= 10) %% Bord Neumann
    invhT=invDiaTri(T1);
    oT=OrdTri(T1);
    CoorBaryT=CoorBary(T1,:);
    [KsipB,Mxy]=KsipDGEdgeB(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF,EdgNormF);
    Ksip(deb1:fin1,deb1:fin1) -= kappa1*KsipB; %% KsipB renvoie -0.5*int_F (...) et pour Neumann il faut 0.5*int_F (...) 
  end
end
%
%------------------------------------------
% MATRICES D'INTEGRATION DANS LES TRIANGLES
%------------------------------------------
%
% MassT(i,j)=     \int_T\phi_i\,phi_j
% KsipT(i,j)=-0.5*\int_T kappa_T* (\Delta\phi_i\,phi_j+\Delta\phi_j\,phi_i)
for t = 1:Nbtri
  oT=OrdTri(t);
  invhT=invDiaTri(t);
  IGLO=NumTri(t,:);
  CoorNeuT=CoorNeu(IGLO,:);
  AGLO=TriEdg(t,:);
  CoorMilT=CoorMil(AGLO,:);
  CoorBaryT=CoorBary(t,:);
  aire=Aires(t);
  %
  MassT=MassDGTriangle(aire,invhT,oT,CoorNeuT,CoorMilT,CoorBaryT);
  debT=idof(t,1); finT=idof(t,2);
  Mass(debT:finT,debT:finT)=MassT;  %
  if (oT>=2)
    KsipT=KsipDGTriangle(invhT,oT,MassT);
    Ksip(debT:finT,debT:finT) += kappa(t).*KsipT;
  end
end


