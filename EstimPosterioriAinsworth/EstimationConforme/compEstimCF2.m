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
% Author : Raphael Lecoq CEA
%
%
% SYNOPSIS [Estim_err_C, Estim_err_C2, norm2_curl_rhok, norm2_rhok, contribution_rhok, norm2_Data, norm2_Diri, norm2_Neu]  = compEstimCF2(idof,Uh)
%   
% Estimateur de l'erreur conforme et quantitées d'intérêt de la p1783 et eq (52)-(53) p1794
% GLOBALS:
%   - kappa(Nbtri,1)         : diffusion dans le triangle 
% INPUT : 
%   - Uh                         : solution numérique 
%   - idof(Nbtri,3)              : index des degrés de liberté
% 
% OUTPUT:
%   - Estim_err_C                : estimation d'erreur conforme
%   - Estim_err_C2(Nbtri,1)      : carré de l'estimation locale de l'erreur conforme
%   - norm2_curl_rhok(Nbtri,1)   : carré de la norme L2 de curl(rho_T)
%   - norm2_rhok(Nbtri,1)        : carré de la norme L2 de rho_T 
%   - contribution_rhok(Nbtri,1) : norm2_curl_rhok + norm2_curl_rhok
%   - norm2_Data(Nbtri,1)        : carré de la norme L2 de l'oscillation des données = norm2_Diri + norm2_Neu
%   - norm2_Neu(Nbtri,1)         : carré de la norme L2 de l'oscillation des données de Neumann
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Estim_err_C, Estim_err_C2, norm2_curl_rhok, norm2_rhok, contribution_rhok2, norm2_Data, norm2_Neu]  = compEstimCF2(idof,Uh)
  
  
global kappa 

[norm2_rhok, norm2_curl_rhok] = norm_rhok(idof,Uh);
contribution_rhok2 = (norm2_rhok -  norm2_curl_rhok)./kappa;

%%% Contibution du second membre 

[norm2_Data,  norm2_Neu] = dataOscillation();
Estim_err_data = sqrt(sum(norm2_Data)  + sum(norm2_Neu));

Estim_err_C2 = (sqrt(norm2_Data) + sqrt(contribution_rhok2)  + sqrt(norm2_Neu)).^2;


Estim_err_C = sqrt(sum(Estim_err_C2));
