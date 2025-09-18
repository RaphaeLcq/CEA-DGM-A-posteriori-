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
function [KsipB,Mxy]=KsipDGEdgeB(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF,EdgNormF)
%
% Matrice de masse arete du bord
Mxy=MassDGEdgeB(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF);
%
%-----------------------
% MATRICE DE CONSISTANCE
%-----------------------
%
ndof=3*oT;
KsipB= zeros(ndof,ndof);
% arete du bord

Svec2  = -invhT*EdgNormF;
Svec   = 0.5*Svec2;
% polynomes d'ordre 0
KsipB(1,2:3)= Svec;
KsipB(2:3,1)= Svec';
% polynomes d'ordre 1
KsipB(2,2)  = Mxy(1,2)*Svec2(1);
KsipB(2,3)  = Mxy(1,3)*Svec(1)+Mxy(1,2)*Svec(2);
KsipB(3,2)  = KsipB(2,3);
KsipB(3,3)  = Mxy(1,3)*Svec2(2);

if (oT>1)
  % polynomes d'ordre 1
  KsipB(1,4)  = KsipB(2,3);
  KsipB(1,5)  = KsipB(2,2);
  KsipB(1,6)  = KsipB(3,3);
  % polynomes d'ordre 2
  KsipB(2,4)  =   Mxy(2,3)*Svec2(1)+Mxy(2,2)*Svec(2);
  KsipB(2,5)  = 3*Mxy(2,2)*Svec(1);
  KsipB(2,6)  =   Mxy(3,3)*Svec(1)+Mxy(2,3)*Svec2(2);
  KsipB(3,4)  =   KsipB(2,6);
  KsipB(3,5)  =   KsipB(2,4);
  KsipB(3,6)  = 3*Mxy(3,3)*Svec(2);
  % bloc symetrique 
  KsipB(4:6,1:3)=KsipB(1:3,4:6)';  
  % polynomes d'ordre 3
  KsipB(4,4)=   Mxy(2,6)*Svec2(1)+  Mxy(2,4)*Svec2(2);
  KsipB(4,5)= 3*Mxy(2,4)*Svec(1) +  Mxy(2,5)*Svec(2);
  KsipB(4,6)=   Mxy(3,6)*Svec(1) +3*Mxy(2,6)*Svec(2);
  KsipB(5,5)= 2*Mxy(2,5)*Svec2(1);
  KsipB(5,6)=   KsipB(4,4);
  KsipB(6,6)= 2*Mxy(3,6)*Svec2(2);
  %
  % symetrie
  KsipB(5,4  )=KsipB(4,5);
  KsipB(6,4:5)=KsipB(4:5,6)';     
end


