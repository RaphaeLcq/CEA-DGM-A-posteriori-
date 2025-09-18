%{
/****************************************************************************
* Copyright (c) 2022, CEA
* All rights reserved.
*
* Redistribution and use iabs(y(i) - 1) <n source and binary forms, with or without modification, are permitted provided that the following conditions are met:
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
% mainPoissonSinusDG.m:
%
% Convergence en maillage, EF DG P1 ou P2, base de monomes ou de Lagrange
%   Pb de Laplace dans un carre [0,1]*[0,1] (Square) ou un Lshape [0;0.5]*[0;1] union [0.5;1]*[0.5;1]
%   -Delta Phi = f
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Paramètres de visu %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global visuEta = 0;    % Permet de différencier le cas où Eta est initialisé avec global EtaEdgMesh ou par la fonction EtaParam.m
global DonneesP1 = 0;  % global DonnesP1 = 1 pour vérifier que \int_T f + \int_F g_T = 0 pour des données de dirichlet projetées sur P1_c.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global DEBUG = 0
global visu = 0 %% 1 cartes locales, 2 solution numérique, 3 si afficher l'erreur L2 aussi
global raffinement = 0 %% Choisi si le maillage est raffiné ou non
global SWIP = 1 %% Choisi si on résout avec SWIP ou SIP 

% Active un problème parmi la liste suivante :
%   'SquareSinus'    : Carré Dirichlet - Δu = nπ² sin(nπx) sin(nπy), u=0 au bord
%   'Neumann'        : Carré Neumann - Δu=0, grad u·n = [3x²-3y² ; 6xy]·n
%   'SquareHarmonic' : À supprimer, mal posé
%   'Lshape'         : L-shape Dirichlet, solution r^α sin(αθ), α=2/3
%   'LshapeNeumann'  : L-shape Neumann, solution r^α sin(αθ), α=2/3
%   'NeumannTop'     : Carré CL mixtes, u=0 sur bottom/left/right, grad u·n=1 sur top
%   'SquareHole'     : Carré avec trou, Dirichlet=0 sauf top, Neumann=-1 top, 0 sur trou
%   'SquareSWIP'     : Problème SquareSWIP
%
problemName =  'SquareSWIP';
setProblem(problemName);


% PARAMETRES MODIFIABLES: mi est le mesh step (1-5) et ef le nombre de simulations (1-2)
m0=1; ef0=1;  
m1=1; ef1=1;
%
global eps nitMAX lambda alpha
eps=1.e-10;  %OUTPUT - KLG, K1, K2 : matrices de raideur locales pour EF de Lagrange ordre, P1 P2 precision GCP
nitMAX=1000; % nb iterations max GCP
lambda = 1;    % Longueur d'onde pour la solution sinusoïdale. Choisir un nombre entier
alpha = 2/3;  % Coefficient alpha pour la solution SWIP
%
nmesh=m1-m0+1; nef=(ef1-ef0)+1; 
ef=[ef0;ef1];
%
meshstep=[0.1,0.05,0.025,0.0125,0.00625];
ndof_tab = zeros(5,1);
TabNbpt  = zeros(1,nmesh);  %[,  568, 2212,  8558,  34239];
TabNbtri = zeros(1,nmesh);  %[242, 1054, 4262, 16794,  68476];
TabNbedg = zeros(1,nmesh);  %[383, 1621, 6473, 25351, 102074];
TabDG    = zeros(nef,nmesh);%[383, 1621, 6473, 25351, 102074];
Er0=zeros(nef,nmesh); Er1=zeros(nef,nmesh);  ErEstim=zeros(nef,nmesh);
Effectivity = zeros(nmesh,1); BrokenNorm=zeros(nef,nmesh);
ErData= zeros(nef,nmesh); Estim_err_C= zeros(nef,nmesh); Estim_err_NC= zeros(nef,nmesh);
Estim_err_C_without_curl = zeros(nef,nmesh); Estim_err_without_curl =  zeros(nef,nmesh);
im=0;
%
% Image finale
%fig=-nmesh*nef;
fig=-nmesh*nef+1;
%
%
global npi npi2
npi=lambda*pi; npi2=2*npi^2;
%
for mii=m0:m1
  im=im+1;
  filename = NomFichierMaillage(mii);
  global mi
  mi=mii;
  %
  global Nbpt CoorNeu CoorNeu2 RefNeu
  global Nbtri NumTri NumTri2 TriEdg CoorBary Aires SomOpp SigTri
  global Nbedg NumEdg CoorMil RefEdg LgEdg2 EdgNorm EdgTri RefTri
  %
  [CoorNeu,CoorNeu2,CoorBary,RefNeu,RefNeu2,NumTri,NumTri2,RefTri,RefTri2,NumEdg,NumEdgB,CoorMil, RefEdg,RefEdgB,TriEdg,EdgTri,SomOpp,LgEdg2,EdgNorm,Aires]=readmeshfiles(filename);
  %
  Nbpt = size(CoorNeu,1); TabNbpt(im) = Nbpt;
  Nbtri= size(NumTri,1) ; TabNbtri(im)= Nbtri;
  Nbedg= size(NumEdg,1) ; TabNbedg(im)= Nbedg;
  %
  global invLgEdg LgEdg
  global DiaTri invDiaTri sigTri
  [DiaTri,invDiaTri,LgEdg,invLgEdg,sigTri]=computeHTri();
  %
  RefTriRegion();
  [D] = paramDiffusion(alpha);
  global kappa
  kappa = kappaTri(D);
  projP1c_Edge(); %% Projection de la donnée de Dirichlet sur P1 continue par morceaux pour l'estimation a posteriori
  %
  ndofDG=[Nbtri,3*Nbtri,6*Nbtri];
  efi=0;
  for i=1:nef
    global ordre
    ordre=ef(i); TabDG(i,im)=ndofDG(ordre+1);
    fig=fig+1;
    efi=efi+1;
    fprintf('Solution P%i.\n',ordre);
    [Er0(efi,im),Er1(efi,im),ErEstim(efi,im),Estim_err_C(efi,im), Estim_err_NC(efi,im),ErData(efi,im),BrokenGradNorm(efi,im),Estim_err_without_curl(efi,im)] = SolvePoissonPourEstimateur();
    Effectivity(mii) = ErEstim(efi,im)/BrokenGradNorm(efi,im);
    ndof_tab(mii) = ndofDG(ordre+1);
  endfor
endfor


if(m0-m1)
%%%%%%%%%%%%%%% ENREGISTREMENT DES ERREURS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_dir = 'resultats_erreurs';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end


savefilename = NomFichierSauvegarde(save_dir, m0, m1);


% Création du tableau des données d'erreur
x_data = sqrt(ndof_tab(m0:m1));

% Ouverture du fichier CSV
fid = fopen(savefilename, 'w');

% Écriture de l'en-tête
fprintf(fid, 'sqrt_ndof,BrokenGradNorm,ErEstim,Effectivity,Estim_err_C,Estim_err_NC,ErData,Estim_err_without_curl\n');

% Écriture des données
for i = 1:length(x_data)
    fprintf(fid, '%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n', ...
        x_data(i), ...
        BrokenGradNorm(1,m0+i-1), ...
        ErEstim(1,m0+i-1), ...
        Effectivity(m0+i-1), ...
        Estim_err_C(1,m0+i-1), ...
        Estim_err_NC(1,m0+i-1), ...
        ErData(1,m0+i-1), ...
        Estim_err_without_curl(1,m0+i-1));
end

fclose(fid);
fprintf('Données d''erreur sauvegardées dans: %s\n', savefilename);


%%%%%%%%%%%%%% VISU DES GRAPHES D'ERREUR

  fig = 50;
  figure(fig);

  x = sqrt(ndof_tab(m0:m1));
  y_left1 = BrokenGradNorm(1,m0:m1);
  y_left2 = ErEstim(1,m0:m1);
  y_right = Effectivity(m0:m1);

  % Axe gauche
  ax1 = axes();
  loglog(ax1, x, y_left1, '-+', 'LineWidth', 2,x, y_left2, '-x','LineWidth', 2);
  set(ax1, 'XScale', 'log', 'YScale', 'log');
  xlabel('sqrt(Ndof) ~ 1/h');
  ylabel(ax1, 'Erreur / Estimation');
  grid(ax1, 'on');
  hold(ax1, 'on');

  % Axe droit superposé
  ax2 = axes('Position', get(ax1, 'Position'), ...
             'XAxisLocation', 'bottom', ...
             'YAxisLocation', 'right', ...
             'Color', 'none', ...
             'XColor', 'none', ...     % <== supprime l'affichage X
             'YColor', 'k', ...
             'XTick', [], ...          % <== supprime les ticks X
             'XLim', get(ax1, 'XLim'));

  % Tracé sur l’axe de droite
  line(x, y_right, 'Color', 'k', 'LineStyle', '-.','Parent', ax2);

  set(ax2, 'XScale', 'log', 'YScale', 'linear');
  ylabel(ax2, 'Effectivité');

  linkaxes([ax1, ax2], 'x');
  title(ax1, "Axe gauche: erreur log-log, Axe droite: effectivité (x-log)");
  hleg = legend(ax1, 'Erreur grad', 'Estimation erreur', 'Location', 'northwest');
  set(hleg, 'FontSize', 14);
1
  hleg2 = legend(ax2, 'Effectivité', 'Location', 'northeast');
  set(hleg2, 'FontSize', 14);


  fig -= 1;
  figure(fig);

  x = sqrt(ndof_tab(m0:m1));
  y1 = ErEstim(1,m0:m1);          % Estimation d’erreur
  y2 = Estim_err_C(1,m0:m1);      % Estimation conforme
  y_ratio = y1 ./ y2;             % Rapport à afficher à droite

  % === Axe principal gauche : deux courbes en log-log ===
  ax1 = axes();
  loglog(ax1, x, y1, '-', 'Color', [0.8500, 0.1250, 0.0500], 'LineWidth', 1.1, ...
               x, y2, '--o', 'MarkerFaceColor', 'red', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.3);

  xlabel('sqrt(Ndof) ~ 1/h');
  ylabel(ax1, "Estimation de l'erreur");
  grid(ax1, 'on');

  % Légende
  hleg = legend(ax1, {'Estimation d’erreur', 'Estimation conforme'}, 'Location', 'northeast');
  set(hleg, 'FontSize', 14);
  title(ax1, "Comparaison des estimateurs d'erreur total et conforme (échelle log-log)");

  % === Auto-zoom sur Y (log scale) ===
  ymin = min([y1(:); y2(:)]);
  ymax = max([y1(:); y2(:)]);
  zoom_margin = 0.1;
  ylim(ax1, [ymin * (1 - zoom_margin), ymax * (1 + zoom_margin)]);

  % === Deuxième axe à droite pour afficher le rapport y1/y2 ===
  ax2 = axes('Position', get(ax1, 'Position'), ...
             'XAxisLocation', 'bottom', ...
             'YAxisLocation', 'right', ...
             'Color', 'none', ...
             'XColor', 'none', ...
             'YColor', [0.2 0.6 0.2], ...
             'XTick', [], ...
             'XLim', get(ax1, 'XLim'));

  % Tracé du rapport (échelle linéaire en Y)
  line(x, y_ratio, 'Color', [0.2 0.6 0.2], 'LineStyle', '--', 'LineWidth', 1.5, 'Parent', ax2);
  set(ax2, 'XScale', 'log', 'YScale', 'linear');
  ylabel(ax2, 'Rapport ErEstim / Estim\_err\_C');

  % === Lien des axes en X ===
  linkaxes([ax1, ax2], 'x');

  % Légende pour le second axe
  hleg2 = legend(ax2, 'Rapport erreur/conforme', 'Location', 'southwest');
  set(hleg2, 'FontSize', 12);

  fig -= 1;
  figure(fig);
  loglog(sqrt(ndof_tab(m0:m1)),BrokenGradNorm(1,m0:m1),'-',sqrt(ndof_tab(m0:m1)),ErEstim(1,m0:m1),'-.',sqrt(ndof_tab(m0:m1)),Estim_err_C(1,m0:m1),'-x', ...
          sqrt(ndof_tab(m0:m1)),Estim_err_NC(1,m0:m1),'-o',sqrt(ndof_tab(m0:m1)),ErData(1,m0:m1),'-+' );
  legend( 'Erreur grad', 'Estimation de l erreur',"Estimation conforme", "Estimation non conforme", "Oscillation données");
  grid on;
  xlabel('sqrt(Ndof) ~ 1/h');
  ylabel('Erreurs');
  title("Erreur et quantités de l'estimation en fonction du ndof, échelle log-log") ;



  fig -= 1;
  figure(fig);
  loglog(sqrt(ndof_tab(m0:m1)),ErEstim(1,m0:m1),'-',sqrt(ndof_tab(m0:m1)),Estim_err_without_curl(1,m0:m1),'-.');
  legend( 'Estimation de l erreur avec la correction du rotationnel', "Estimation de l'erreur sans la correction du rotationnel");
  grid on;
  xlabel('sqrt(Ndof) ~ 1/h');
  ylabel('Erreurs');

  if visu == 3
    fig -= 1;
    figure(fig);
    loglog(sqrt(ndof_tab(m0:m1)),Er0(1,m0:m1));
    legend("Erreur L2 SIP");
    grid on;
    xlabel('sqrt(Ndof) ~ 1/h');
    ylabel('Erreurs');
  endif
endif
