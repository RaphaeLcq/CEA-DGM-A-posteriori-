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
% DGXYtoLG.m: 
%      
% INPUT  - Phih(Ndof,1) : vecteur DG base des monomes
%%       - ordreLG : ordre de la projection
%        - idof(Nbtri,3)  : numero globaux
%
% OUTPUT - ULG(ndofLG,1) : valeurs du vecteur Phih aux sommets (ordreLG=1) 
%                           ou aux sommets et aux aretes (ordreLG=2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ULG=DGXYtoLG(Phih,ordreLG,OrdTri,MLG,idof)
%
global Nbedg CoorMil TriEdg
global Nbtri NumTri TriEdg CoorBary Aires invDiaTri
global Nbpt CoorNeu
%
NdofLG=3*ordreLG;
Nlg=Nbpt;
if (ordreLG==2)
  Nlg=Nlg+Nbedg;
end
%
PUh=zeros(Nlg,1);
%
for t=1:Nbtri
  invhT=invDiaTri(t);
  aire=Aires(t);
  NdofDG=idof(t,3);
  ordreDG=OrdTri(t,1);
  %
  % matrice de projection monone vers Pk
  % PUh(i)=int_T \phi_i\,\psi_j
  % avec \phi_i : base Pk et \psi_j base monomiale
  %
  Mloc=zeros(NdofLG,NdofDG);
  IGLO=NumTri(t,:);
  CoorNeuT=CoorNeu(IGLO,:);
  [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
  CoorBaryT=ones(np,1)*CoorBary(t,:);
  xypT=invhT*(xyp-CoorBaryT);
  awp=aire*wp';
  phiDG=[ones(np,1),xypT];
  if (ordreDG==2)
    phiDG=[phiDG,xypT(:,1).*xypT(:,2),xypT.*xypT];
  end
  debT=idof(t,1); finT=idof(t,2);
  PhihT=awp.*(phiDG*Phih(debT:finT,1));
  if (ordreLG==1)
    psiLG=lambda;
  else
    AGLO=TriEdg(t,:); 
    JLOC=[2,3,1]; KLOC=[3,1,2];
    lambda2=lambda.*lambda;
    psiS=2*lambda2-lambda;
    psiM=4*lambda(:,JLOC).*lambda(:,KLOC);
    psiLG=[psiS,psiM]; IGLO=[IGLO,AGLO+Nbpt]; 
  end 
  PUh(IGLO)+=sum(PhihT.*psiLG,1)';
end % end_for t=1:Nbtri
ULG=MLG\PUh;