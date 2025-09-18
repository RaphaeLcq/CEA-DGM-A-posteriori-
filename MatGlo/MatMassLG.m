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
% MatMassLG.m:
% Assemblage matrice de masse EF de Lagrange P1 ou P2
%
% SYNOPSIS MLG = MatMassLG(ordre)
%          
% INPUT  - Nbpt : nombre de sommets
%        - Nbedg : nombre d'aretes
%        - NumTri(Nbtri,3) : liste de triangles 
%                   (3 numeros de sommets)
%        - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'Edg opposee au sommet NumTri(l,i)
%                  (3 numeros des Edg - matrice entiere Nbtri x 3)
%		     - Aires(Nbtri,1) : aires des triangles
%        - ordre : 0, 1 ou 2
%
% OUTPUT - MLG : matrice de masse 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MLG = MatMassLG(ordre)
%
global Nbtri NumTri TriEdg Aires
global Nbpt
global Nbedg
%
% Matrices de masse locales pour aire=1 !!
[MT1,invMT1,MT01,MT2,invMT2,MT02,MT12]=MassLGTriangle();
%
% ORDRE 0
if (ordre==0)
  MLG=diag(sparse(Aires),0);
else
  % ORDRE 1
  if (ordre==1)
     MT=MT1; Ndof=Nbpt; nloc=3; Numloc=NumTri;
  end
  % ORDRE 2
  if (ordre==2)
     MT=MT2; Ndof=Nbedg+Nbpt; nloc=6; Numloc=[NumTri,Nbpt+TriEdg];
  end
  MLG=sparse(Ndof,Ndof);
  for t=1:Nbtri
     MLG_loc=Aires(t)*MT;
     IGLO=Numloc(t,:);
     MLG(IGLO,IGLO)+=MLG_loc;
  end % end_for t=1:Nbtri
end
