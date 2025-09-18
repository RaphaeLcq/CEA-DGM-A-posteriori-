%{
/****************************************************************************
* Copyright (c) 2023, CEA
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
% GCPDG.m:
%
% Algorithme de resolution du systeme mixte 
% A U=F
%
% SYNOPSIS Uh=GCPDG(Ku,RHS_U)
%          
% Uzawa : on resout A U=F par un GCP par blocs
%
% INPUT : - Ku(ndofU,ndofU)  : matrice de raideur interne (une composante)
%         - RHS_U(ndofU,1)     : second membre 
% OUTPUT: - Uh(ndofU,1) : solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Uh,res,nit]=GCPDG(invKu,Ku,RHS_U)
%
global eps nitMAX
%
%-------------------------%
% Initialisation AUh=Fu GCP
%-------------------------%
Uh=invKu*RHS_U;
%----%
% GCP
%----%
r=RHS_U-Ku*Uh;
z=invKu*r;
p=z;
rr=r'*r;
eps2=eps^2*rr;
rz=r'*z;
nit=0;
while ((rr>eps2)&&(nit<nitMAX))
  Ap=Ku*p;
  alfa=(r'*z)/(Ap'*p);
  Uh+=alfa*p;
  r-=alfa*Ap;
  z=invKu*r;
  rznew=(r'*z);
  beta=rznew/rz;
  rz=rznew;
  p=z+beta*p;
  rr=r'*r;
  nit+=1;
end
res=r;