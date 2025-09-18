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
%
% Matrices de raideur et matrice de masse
% Formulation DG SIP, fonctions de base monomiales
%
% SYNOPSIS [Ksip,Mass]=KsipPolyDGsip(OrdTri,idof)
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
%
% OUTPUT - Ksip : matrice  de raideur, associee a l'operateur -Delta u, fonctions de base monomiales
%        - Mass : matrice de masse, fonctions de base monomiales
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [Ksip,Mass]=KsipPolyDGsip(OrdTri,idof)
%
global CoorNeu
global Nbtri CoorBary Aires NumTri TriEdg invDiaTri
global Nbedg NumEdg CoorMil RefEdg EdgNorm EdgTri
global etaEdg CoorBary LgEdg EdgC kappa SWIP

%
Ndof=idof(Nbtri,2);
Ksip=sparse(Ndof,Ndof); 
Mass=sparse(Ndof,Ndof);

%
%---------------------------------
% MATRICES D'INTEGRATION AUX FACES
%---------------------------------
%;
for f=1:Nbedg
  ij=NumEdg(f,:);
  CoorEdgF=CoorNeu(ij,:);
  CoorMilF=CoorMil(f,:);
  EdgNormF=EdgNorm(f,:);
  T12=EdgTri(f,:);
  T1=T12(1); deb1=idof(T1,1); fin1=idof(T1,2); 
  etaF=etaEdg(f);
  LF=LgEdg(f,1);
  etaK = etaKappa(f);
  switch SWIP
    case 0
      [kappaAvg1, kappaAvg2] = [0.5,0.5];
    case 1
      [kappaAvg1, kappaAvg2] = HarmonicAverageKappa(f);
  endswitch

  if (RefEdg(f)==0)
    oT=OrdTri(T12);
    T2=T12(2); deb2=idof(T2,1); fin2=idof(T2,2);
   if ((oT(1)==0)&&(oT(2)==0))              
     Ksip(deb1:fin1,deb2:fin2)=-1;                
     else     
     
    % Face interieure
    % Kxy12(i,j)= -0.5*\int_F\grad\psi_j\cdot\nvec_F\,\phi_i  CONSISTANCE
    %             +0.5*\int_F\grad\phi_i\cdot\nvec_F\,\psi_j  SYMETRIE
    % Mxy1(i,j) =  \int_F \phi_i \phi_j
    % Mxy2(i,j) =  \int_F \psi_i \psi_j
    % Mxy12(i,j)=  \int_F \phi_i \psi_j
    invhT=invDiaTri(T12);
    oT=OrdTri(T12);
    CoorBaryT=CoorBary(T12,:);
    IGLO1=NumTri(T1,:);
    CoorNeuT1=CoorNeu(IGLO1,:);
    IGLO2=NumTri(T2,:);
    CoorNeuT2=CoorNeu(IGLO2,:);
    kappa1 = kappa(T1); 
    kappa2 = kappa(T2);
    [Kxy1,Kxy2,Kxy12,SVxy1,SVxy2,SVxy12]=KsipDGEdgePoly(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF,EdgNormF,LF,kappaAvg1,kappaAvg2,kappa1,kappa2);
    Ksip(deb1:fin1,deb1:fin1)+=-Kxy1+etaF*etaK*SVxy1;
    Ksip(deb2:fin2,deb2:fin2)+=-Kxy2+etaF*etaK*SVxy2;
    Ksip(deb2:fin2,deb1:fin1)+=-Kxy12+etaF*etaK*SVxy12;
  end
        % SYMETRIE
    Ksip(deb1:fin1,deb2:fin2)=Ksip(deb2:fin2,deb1:fin1)';
  else
    % Face du bord
    % KsipB(i,j)=-0.5*\int_F\grad(\phi_i\,\phi_j)\cdot\nvec_F
    %            +etaF\int_F\phi_i\,\phi_j
    %Ksip(deb1:fin1,deb1:fin1)=Ksip(deb1:fin1,deb1:fin1)-1; 
    invhT=invDiaTri(T1);
    oT=OrdTri(T1);
    CoorBaryT=CoorBary(T1,:);
    IGLO=NumTri(T1,:);
    CoorNeuT=CoorNeu(IGLO,:);
    if (oT(1)~=0)
    [KsipB,SVxy]=KsipDGEdgeBPoly(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF,EdgNormF,LF);
    if  (RefEdg(f)>0 & RefEdg(f) < 10)
       Ksip(deb1:fin1,deb1:fin1)+= -kappa(T1)*KsipB+kappa(T1)*etaF*SVxy;
     end
   end
  end
end
%
%------------------------------------------
% MATRICES D'INTEGRATION DANS LES TRIANGLES
%------------------------------------------
%
% MassT(i,j) = \int_T\phi_i\,phi_j
% KsipT(i,j) = \grad phi_i \grad phi_j
for t=1:Nbtri
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

  KT=RaidPolyTriangle(aire,invhT,oT,CoorNeuT,CoorMilT,CoorBaryT);
  debT=idof(t,1); finT=idof(t,2);
  Ksip(debT:finT,debT:finT)+=kappa(t)*KT;
  %
  if oT==0
    Ksip(debT:finT,debT:finT)=kappa(t)*3;
    end
end



