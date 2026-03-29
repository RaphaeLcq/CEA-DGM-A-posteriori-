
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
% Author : CEA
%
%
% SYNOPSIS
% Computes the weighted average on each edge.
%
% GLOBALS:
%   - Nbtri  : global number of triangles
%   - RefEdg(Nbtri,1) : reference of edges of the mesh
%   - kappa(Nbtri,1)   : diffusion in the triangle
%   - SWIP : global of SWIP method, 1 if using SWIP, 0 if using SIP
% INPUT :
%   - edge : selected edge
%
% OUTPUT:
%   - avg1, avg2   : harmonic averages on each triangle sharing the edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [avg1,avg2] = HarmonicAverageKappa(edge)


  global EdgTri RefEdg kappa SWIP
  avg1 = 0; avg2 = 0;
  if (SWIP==1)
    if (RefEdg(edge) > 0) %% Dirichlet boundary
      tri = EdgTri(edge,1);
      avg1 = kappa(tri);
      avg2 = 0;
    else
      tri1 = EdgTri(edge,1); tri2 = EdgTri(edge,2);
      kappa1 = kappa(tri1); kappa2 = kappa(tri2);
      avg1 = (kappa2)/(kappa1 + kappa2);
      avg2 = (kappa1)/(kappa1 + kappa2);
    endif
  else
    avg1 = 0.5;
    avg2 = 0.5;
  endif

endfunction
