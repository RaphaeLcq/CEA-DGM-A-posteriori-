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
% ValPolyEdge.m:
% 
% DG-FEM (SIPG monomes P1 ou P2)
%
% SYNOPSIS [wp3,Mat3T1,Mat3T2,wp5,Mat5T1,Mat5T2]=ValPolyEdge(invhT,oT,CoorBaryT,CoorEdgF)
%          
% INPUT : - invhT  : mise a l'echelle locale (inverse du diametre du triangle)
%         - oT(2,1)   : ordres d'approximation locaux aux triangles
%         - CoorBaryT(1,2): coordonnees des barycentres des triangles
%         - CoorEdgF(2,2) : coordonnees de l'arete
% OUTPUT: - wp3   : poids integration sur arete a 3 points
%         - Mat3T1: valeurs polynomes pour integration a 3 points T1
%         - Mat3T2: valeurs polynomes pour integration a 3 points T2
%         - wp5   : poids integration sur arete a 5 points
%         - Mat5T1: valeurs polynomes pour integration a 5 points T1
%         - Mat5T2: valeurs polynomes pour integration a 5 points T2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wp3,Mat3T1,Mat3T2,wp5,Mat5T1,Mat5T2]=ValPolyEdge(invhT,oT,CoorBaryT,CoorEdgF)
  %
  [xyp3F,wp3,lambda3,np3]=IntEdg_Sim3(CoorEdgF);
  xyp3T1=invhT(1)*(xyp3F-ones(np3,1)*CoorBaryT(1,:));
  Mat3T1=[xyp3T1,xyp3T1(:,1).*xyp3T1(:,2),xyp3T1.^2];
  %
  xyp3T2=invhT(2)*(xyp3F-ones(np3,1)*CoorBaryT(2,:));
  Mat3T2=[xyp3T2,xyp3T2(:,1).*xyp3T2(:,2),xyp3T2.^2];
  %
  o1o2=oT(1)+oT(2);
  wp5=[]; Mat5T1=[]; Mat5T2=[];
  if (o1o2>2)
    [xyp5F,wp5,lambda5,np5]=IntEdg_Boo5(CoorEdgF);
    if (oT(1)>1)
      xyp5T1=invhT(1)*(xyp5F-ones(np5,1)*CoorBaryT(1,:));
      Mat5T1=[xyp5T1,xyp5T1(:,1).*xyp5T1(:,2),xyp5T1.^2];
    end
    %
    if (oT(2)>1)
      xyp5T2=invhT(2)*(xyp5F-ones(np5,1)*CoorBaryT(2,:));
      Mat5T2=[xyp5T2,xyp5T2(:,1).*xyp5T2(:,2),xyp5T2.^2];
    end
  end
end % end_function