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
% Author :  Erell Jamelot CEA
%
% ValPolyEdgeB.m:
% 
% DG-FEM (SIPG monomes P1 ou P2)
% 
% Valeur des polynomes aux points d'integration sur une arete
%
% SYNOPSIS [wp3,Mat3T,wp5,Mat5T]=ValPolyEdgeB(invh,ordre,CoorBary,CoorEdg)
%          
% INPUT : - invhT  : mise a l'echelle locale (inverse du diametre du triangle)
%         - oT : ordre d'approximation local au triangle
%         - CoorBaryT(1,2): coordonnees du barycentre du triangle
%         - CoorEdgF(1,2) : coordonnees de l'arete
% OUTPUT: - wp3  : poids integration sur arete a 5 points
%         - Mat3T: valeurs polynomes pour integration a 3 points
%         - wp5  : poids integration sur arete a 5 points
%         - Mat5T: valeurs polynomes pour integration a 5 points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wp3,Mat3T,wp5,Mat5T]=ValPolyEdgeB(invhT,oT,CoorBaryT,CoorEdgF)
  %
  % Schemas d'integration
  [xyp3F,wp3,lambda3,np3]=IntEdg_Sim3(CoorEdgF);
  xyp3T=invhT*(xyp3F-ones(np3,1)*CoorBaryT);
  Mat3T=[xyp3T,xyp3T(:,1).*xyp3T(:,2),xyp3T.*xyp3T];
  %
  wp5=[]; Mat5T=[];
  if (oT>1)
    [xyp5F,wp5,lambda5,np5]=IntEdg_Boo5(CoorEdgF);
    xyp5T=invhT*(xyp5F-ones(np5,1)*CoorBaryT);
    Mat5T=[xyp5T,xyp5T(:,1).*xyp5T(:,2),xyp5T.*xyp5T];
  end
end % end_function