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
% SolutionPoissonSinusLG.m:
% Solution du probleme de Poisson "Sinus"
%
%   Pb de Poisson dans un carre [0,1]*[0,1] :
%   -Delta Phi =(2*npi^2)*sin(npi*x)*sin(npi*y)
%   Solution Phi(x,y)=sin(npi*x)*sin(npi*y)
%
%
% SYNOPSIS [Phiex]=SolutionPoissonSinusLG
%
% GLOBAL - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets (noeuds P1)
%        - CoorBary(Nbtri,2)   : Coordonnees des barycentres de triangles
%        - invDiaTri(Nbtri,2)   : inverse des diametres des triangles
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets) 
%		     - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%        - Aires(Nbtri,1) : aires des triangles
% INPUT  - ordre  : ordre d'approximation
%        - MLG : matrice de masse
% OUTPUT - [PhiexLG] : la solution exacte au points de discretisation
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PhiexLG]=SolutionPoissonSinusLG(ordre, MLG)
%
global Nbpt CoorNeu 
global Nbtri NumTri Aires TriEdg
global Nbedg
global npi
%
Ndof=Nbpt;
if (ordre==2)
  Ndof+=Nbedg;
end
RHS=zeros(Ndof,1); Phiex=zeros(Ndof,1);
for t=1:Nbtri
  IGLO=NumTri(t,:); CoorNeuT=CoorNeu(IGLO,:);
  %   
  % volume of the element
  aire=Aires(t);
  % Integration points
  [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
  awp=aire*wp';
  sinxy=sin(npi*xyp);
  PhiexT=awp.*sinxy(:,1).*sinxy(:,2);
  if (ordre==1)
    psiLG=lambda;  UGLO=IGLO;
  else
    AGLO=TriEdg(t,:); 
    JLOC=[2,3,1]; KLOC=[3,1,2];
    lambda2=lambda.*lambda;
    psiS=2*lambda2-lambda;
    psiM=4*lambda(:,JLOC).*lambda(:,KLOC);
    psiLG=[psiS,psiM]; UGLO=[IGLO,AGLO+Nbpt]; 
  end
  %
  RHS(UGLO,1)+=sum(PhiexT.*psiLG,1)';
end
%
PhiexLG=MLG\RHS;
