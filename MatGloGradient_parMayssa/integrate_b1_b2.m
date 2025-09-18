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
% Author : Mayssa Mroueh CEA
%
%
% SYNOPSIS
    % Calcule l'intégrale de (b1.(x-x1)) * (b2.(x-x2)) sur un segment en 2D.
    %
    % Entrées :
    % - b1 : vecteur directionnel pour le premier terme [bx1, by1]
    % - x1 : point associé à b1 [x1, y1]
    % - b2 : vecteur directionnel pour le second terme [bx2, by2]
    % - x2 : point associé à b2 [x2, y2]
    % - P0 : point de départ du segment [x0, y0]
    % - P1 : point d'arrivée du segment [x1, y1]
    %
    % Sortie :
    % - I : valeur exacte de l'intégrale.

    % Calcul de la longueur et des vecteurs directionnels



function I = integrate_b1_b2(b1, x1, b2, x2, P0, P1)
    L = norm(P1 - P0);
    dP = P1 - P0;

    % Calcul des coefficients
    A = dot(b1, P0 - x1) * dot(b2, P0 - x2);
    B = dot(b1, P0 - x1) * dot(b2, dP) + dot(b2, P0 - x2) * dot(b1, dP);
    C = dot(b1, dP) * dot(b2, dP);

    % Intégrale
    I = L * (A + B / 2 + C / 3);
end