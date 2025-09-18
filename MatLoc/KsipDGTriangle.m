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
% KsipDGTriangle.m:
%
% DG-FEM (SIPG monomes P1 ou P2)
% 
% Matrices de raideur locales, fonctions de base monomiales pour un triangle donne
%
% $\phi_i$ base locale dans $T$ (test fonction)
% $\phi_j$ base locale dans $T$ (trial function),
% $\vec{n}_F$ oriente vers l'exterieur de $T$
%
% Ksip(i,j)= -0.5*(\int_T\Delta\phi_i\,\phi_j+\Delta\phi_j\,\phi_i)
%          = \int_T\grad\phi_i\cdot\grad\phi_j-0.5\sum_{F\in\pa T}\int_F\grad(\phi_i\,\phi_j)\cdot\nvec
%
% 
% SYNOPSIS Ksip=KsipDGTriangle(invh,ordre,Mass)
%          
% INPUT  - invhT  : mise a l'echelle (inverse du diametre du triangle)
%        - oT : ordre d'approximation
%        - MassT  : matrice de masse locale
% OUTPUT - Ksip : matrice associee a l'operateur -Delta u, fonctions de base monomiales
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ksip=KsipDGTriangle(invhT,oT,MassT)

ndof=3*oT;
Ksip=zeros(ndof,ndof);
%
%-----------------------
% MATRICE DE CONSISTANCE
%-----------------------
%
coeff=-(invhT^2);
if (oT>1)
  Ksip(:,5) =  coeff*MassT(:,1);
  Ksip(:,6) =  Ksip(:,5);
  Ksip(5,6) =  Ksip(5,6)+Ksip(6,6);
  Ksip(5,5) = 2*Ksip(5,5);
  Ksip(6,6) = 2*Ksip(6,6);
  % SYMETRIE
  Ksip(5:6,1:4) = Ksip(1:4,5:6)';
  Ksip(6,5)     = Ksip(5,6);
end