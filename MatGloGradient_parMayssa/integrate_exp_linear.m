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
    % Calcule l'intégrale de exp(b1.(x - x1)) * (b3.(x - x3)) sur un segment en 2D.
    %
    % Entrées :
    % - b1 : vecteur directionnel pour l'exponentielle [bx1, by1]
    % - x_1 : point associé à b1 [x*, y*]
    % - b3 : vecteur pour le produit scalaire linéaire [bx3, by3]
    % - x3 : point associé à b3 [x3, y3]
    % - P0 : point de départ du segment [x0, y0]
    % - P1 : point d'arrivée du segment [x1, y1]
    %
    % Sortie :
    % - I2 : valeur exacte de l'intégrale.

    % Calcul de la longueur et des vecteurs directionnels



function I2 = integrate_exp_linear(b1, x_1, b3, x3, P0, P1)
    L = norm(P1 - P0);
    dP = P1 - P0;

    % Coefficients pour l'exponentielle
    c0 = dot(b1, P0) - dot(b1, x_1);
    c1 = dot(b1, dP);
    % Points pour le produit scalaire linéaire
    linear_coeff = dot(b3, P0 - x3);
    linear_grad = dot(b3, dP);
    % Intégrale analytique
    if abs(c1) > 1e-12
        exp_integral = L*(exp(c0+c1)-exp(c0))/c1;
        I2 = exp_integral * linear_coeff + ...
             L*(linear_grad*((c1-1)*exp(c1+c0)+exp(c0))/c1^2);
    else
        exp_integral = L*exp(c0);
        I2 = exp_integral*(linear_coeff+linear_grad/2);
    end
%    [xp,wp,lambdaedg,npedg,rhoedg,thetaedg]=IntEdg_Boo5([P0;P1]);
%    wp=L*wp;
%    I2-wp*(exp((xp-x_1)*b1').*((xp-x3)*b3'))
end
