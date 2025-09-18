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
%%
% SYNOPSIS etaEdg=EtaParam(ordrU)
%          Calcul du parametre eta from Lemma 4.12 p. 129 DiEr12 et WaHe03
%          Pour l'operateur Laplacien $\Delta$ (sinon, tenir compte de nu)
% GLOBAL - Nbtri  : nombre de triangle
%        - Aires(Nbtri,1) : aires des triangles
%        - sigTri(Nbtri,1)    : sigma du triangle (je sais pas ce que c'est, on l'utilise dans EtaParam.m, peut être le nombre de faces ?) 
%        - Nbedg : nombre d'aretes
%        - LgEdg2(Nbedg,1) : longueurs des aretes au carre
%		     - RefEdg(Nbedg,1) : Reference de chaque arete 
%		     - EdgTri(Nbedg,2) : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere

% INPUT  - ordrU(Nbtri,1)  : ordre d'approximation dans chaque triangle
% OUPUT  - etaEdg = (eta/h_F)*|F| parametre de stabilite (cas ou nu=1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function etaEdg=EtaParam(ordrU)
  global Nbtri Aires sigTri
  global Nbedg LgEdg2 RefEdg EdgTri
  %global nu
  %
  etaEdg=zeros(Nbedg,1);
##  etaTri=3*(2/pi)*(sigTri.^2).*(ordrU+1).*(ordrU+2);
  etaTri=(2/pi)*(sigTri.^2).*(ordrU+1).*(ordrU+2);
  % (ordrU+1).*(ordrU+2)=6*ordrU % pour ordrU=1,2
##  etaTri=9*ordrU.*(4/pi)*(sigTri.^2);
  for e=1:Nbedg
    T12=EdgTri(e,:);
    etaT1=etaTri(T12(1));
    if (RefEdg(e)==0)
      etaT2=etaTri(T12(2));
      etaEdg(e)=2*etaT1*etaT2/(etaT1+etaT2);
    else
      etaEdg(e)=etaT1;
  end
  end

endfunction