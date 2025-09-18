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
% Author : Raphaël Lecoq CEA
%
%
% SYNOPSIS - Evalue fonction de Neumann du problème considéré en (x,y)
%

% GLOBAL - Lshape SquareHarmonic Neumann NeumannTop LshapeNeumann SquareHole SquareSinus SquareSWIP  : problème considéré
% Inputs - x(np,1) 
%        - y(np,1) 
%        - edgeGlo : Indice global de l'arête
% Output - value(np,1) : valeur de la donnée de Neumann en ce point


function value = eval_donneesNeumann(x,y,edgeGlo)
  
  global Lshape SquareHarmonic Neumann NeumannTop LshapeNeumann SquareHole SquareSinus SquareSWIP
  
  global CoorNeu NumEdg EdgUnit
  np = length(x);
  value = zeros(np,1);
  if (SquareSinus == 1)
    value = zeros(np,1); 
    
  elseif (Lshape == 1)
    value = zeros(np,1); 
    
  elseif (SquareHarmonic == 1)
    value = zeros(np,1); 
    
  elseif (Neumann == 1)
    global npi
    CN = CoorNeu(NumEdg(edgeGlo,:),:);
    N = EdgUnit(edgeGlo, :);  % N = [nx, ny]
    dfdx = npi.*cos(npi * x).*sin(npi .*y);
    dfdy = npi.*sin(npi * x).*cos(npi .*y);
    value = dfdx .* N(1) + dfdy .* N(2);
    
  elseif (NeumannTop == 1)
    CN = CoorNeu(NumEdg(edgeGlo,:),:);
    if abs(CN(1,2) - CN(2,2)) == 0 %% edge horizontal x == cste
      for i = 1:np
        if (abs(y(i) - 1) < 1e-10)
          value(i) = 1;
        endif
      endfor
    endif
    
  elseif LshapeNeumann == 1
    global LgEdg
    unTiers = 1/3;
    CN = CoorNeu(NumEdg(edgeGlo,:),:);
    N = EdgUnit(edgeGlo, :);  % N = [nx, ny]
    [rp, thetap] = CartesianToPolarCentered(x, y);
    df_dx = 2 * unTiers  .* rp.^(-unTiers) .* (sin(2 * unTiers .* thetap).*cos(thetap) - cos(2 * unTiers * thetap) .*sin(thetap));
    df_dy = 2 * unTiers  .* rp.^(-unTiers) .* ( cos(2 * unTiers * thetap).*cos(thetap) + sin(2 .* unTiers * thetap ).*sin(thetap));
    
    % Vérifier si le point (1/2, 1/2) existe et mettre grad f = avg(grad f, boule de rayon 0.1*h) pour cet indice
    tolerance = 1e-12;  % Tolérance pour la comparaison des flottants
    idx_special = find(abs(x - 1/2) < tolerance & abs(y - 1/2) < tolerance);
    if ~isempty(idx_special)
        df_dx(idx_special) = -(0.1*LgEdg(edgeGlo))^(-1/3)*9/(5*pi);
        df_dy(idx_special) = (0.1*LgEdg(edgeGlo))^(-1/3)*9*sqrt(3)/(10*pi);
    end
    value = df_dx .* N(1) + df_dy .* N(2);
    
  elseif (SquareHole == 1) %% -1 au top, 0 sur le trou carré du milieu
    CN = CoorNeu(NumEdg(edgeGlo,:),:);
    if abs(CN(1,2) - CN(2,2)) == 0 %% edge horizontal x == cste
      for i = 1:np
        if (abs(y(i) - 1) < 1e-10)
          value(i) = 1;
        endif
      endfor
    endif
  
  elseif (SquareSWIP == 1)
    value = zeros(np,1);
  endif
      
  
endfunction
