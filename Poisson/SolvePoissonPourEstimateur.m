%{
/****************************************************************************
* Copyright (c) 2025, CEA
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
% Author : Raphael Lecoq, Erell Jamelot CEA
%
% SolvePoissonPourEstimateur.m:
%
% Mesh convergence, Finite Elemenst Discontinuous Galerkin P1 or P2, basis X^n or Lagrange
%
% SYNOPSIS [Eu0,Eu1,fig]=SolvePoissonPourEstimateur(fig)
%
% GLOBAL
%        - CoorNeu(Nbpt,2)        : coordinates (x, y) of vertices (nodes P1)
%        - CoorNeu2(Nbpt+Nbedg,2) : coordinates (x, y) of nodes P2
%        - RefNeu(Nbpt,1)         : reference of vertices
%        - Nbtri                  : global number of triangles
%        - CoorBary(Nbtri,2)      : coordinates of barycentres of triangles
%        - Aires(Nbtri,1)         : area of triangles
%        - NumTri(Nbtri,3)        : list of triangles  (3 vertices number)
%        - NumTri2(4*Nbtri,3)     : list of triangles of the mesh P2 (3 number of vertices)
%		     - TriEdg(Nbtri,3)        : for each triangle, TriEdg(l,i) is the number of vertice opposite from vertice NumTri(l,i)
%                                  (3 number of edges - full matrix Nbtri x 3)
%        - invDiaTri(Nbtri,2)     : inverse of diametres of triangles
%        - Nbedg                  :  global number of edges
%        - CoorMil(Nbedg,2)       : coordinates of the edges center
%		     - RefEdg(Nbedg,1)        : reference of each vertice
%        - LgEdg2(Nbedg,1)        : length of edges squared
%        - EdgNorm(Nbedg,2)       : vector normal-edge, oriented tri1->tri2
%		     - EdgTri(Nbedg,2)        : for each vertice, EdgTri(a,:) gives the number of the 2 triangles sharing vertice
%                                 EdgTri(a,2) = 0 if a is on the domain border
%        - ordre                  :  approximation order
%        - mi                     : number of the mesh
%        - etaEdg(NbEdg,1)        : value of eta on each edge
% OUTPUT
%        - Er0                    : L2 Error  \|\phi_h-\Pi_h(\phi)\|_0
%        - Er1                    : H1 Error    |\phi_h-\Pi_h(\phi))|_1
%        - ErEstim                : Estimation for Broken Gradient norm
%        - Estim_err_C            : Estimation of the conform error
%        - Estim_err_NC           : Estimation of non conformity (NC) error
%        - Estim_err_data         : Data oscillation
%        - GradNorm               : Broken gradient error sqrt( \sum\limits_T | grad U_ex - grad U |_L² )
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Er0,Er1,Estim_err,Estim_err_C, Estim_err_NC,Estim_err_data,GradNorm, Estim_err_without_curl] = SolvePoissonPourEstimateur()
%
global CoorNeu CoorNeu2 RefNeu
global Nbtri CoorBary Aires NumTri NumTri2 TriEdg invDiaTri
global Nbedg CoorMil RefEdg LgEdg2 EdgNorm EdgTri
global ordre mi
global etaEdg
global visu
%
OrdTri = ordre * ones(Nbtri,1);

% Number of degrees of liberty
ndof = 3 * OrdTri;
idof = zeros(Nbtri,3);
idof(:,3) = ndof;
ndof_tot = sum(ndof);
idof(1,1) = 1;
idof(1,2) = ndof(1);
for t = 2:Nbtri
  idof(t,1) = idof(t-1,2) + 1;
  idof(t,2) = idof(t-1,2) + ndof(t);
end

% Renumerotation
Profil = sparse(Nbtri,Nbtri);
for e = 1:Nbedg
  tri = EdgTri(e,:);
  if (tri(2) > 0)
    Profil(tri,tri) += ones(2,2);
  end
end
renumT = symrcm(Profil);
sU = zeros(ndof_tot,1);
ideb = 1;
for tt = 1:Nbtri
  t = renumT(tt);
  i1 = idof(t,1); ndoft = idof(t,3);
  for p = 0:ndoft-1
    sU(ideb+p,1) = i1 + p;
  end
  ideb += ndoft;
end

% Stability parameter Eta. If VisuEta global is 1, it allows to display the value of the error and estimation for different values of Eta.
global visuEta
if (visuEta)
  global EtaEdgMesh
  etaEdg = EtaEdgMesh;
else
  etaEdg = EtaParam(OrdTri);
end

% MATRICES
fprintf(' Construction of matrices.\n');
[Ku,Mu] = MatPoissonDGsip(OrdTri,idof);

invMu = sparse(ndof_tot,ndof_tot);
for t = 1:Nbtri
  debT = idof(t,1); finT = idof(t,2);
  invMu(debT:finT,debT:finT) = inv(Mu(debT:finT,debT:finT));
end

%  Right hand side and exact solution
fprintf(' Construction of right hand side.\n');
global Lshape SquareHarmonic DonneesP1
if (DonneesP1 == 1 && (Lshape + SquareHarmonic == 1))
  [Phiexh,NormPhiex2h,NormGPhiex2h,RHS] = SolutionPoissonSinusDGLshape_DonneesP1(Mu,invMu,Ku,idof);
else
  [Phiexh,NormPhiex2h,NormGPhiex2h,RHS] = SolutionPoissonDG(Mu,invMu,Ku,idof);
end

% Check the system found the right solution
VERIF = 1;
if (VERIF)
  fprintf(' Checking linear system sanity.\n');
  errNum = abs((Ku*Phiexh - RHS)' * Phiexh) / NormGPhiex2h;
  fprintf('P%i DG mesh_%i, errNum = %7.2e\n', ordre, mi, errNum);
  fprintf('------------------------------------------------\n');
end

% Solving
fprintf(' Solving.\n');
id = tic;
Ku_s = Ku(sU,sU);
RHS_s = RHS(sU,1);

is_neumann_pur = isempty(find(RefEdg > 0 & RefEdg < 10));
if is_neumann_pur
  one_vec = ones(ndof_tot,1);
  m1_s = Mu(sU,:) * one_vec;  % m1_s = Mu_s * 1
  Ku_s_aug = [Ku_s, m1_s; m1_s', 0];
  RHS_s_aug = [RHS_s; 0];
else
  Ku_s_aug = Ku_s;
  RHS_s_aug = RHS_s;
end

% Direction solver with no preconditioning
GCP = 0; CH = 0;
if (GCP == 1)
  invKu_s = sparse(ndof_tot, ndof_tot);
  for t = 1:Nbtri
    i1 = idof(t,1); i2 = idof(t,2);
    i1_s = sU(i1); i2_s = sU(i2);
    Kuloc = Ku(i1:i2, i1:i2);
    invKu_s(i1_s:i2_s, i1_s:i2_s) = inv(Kuloc);
  end
  [Phih_s, res, nit] = GCPDG(invKu_s, Ku_s_aug, RHS_s_aug);
  fprintf('Nb it GCP = %i.\n', nit);
else
  if (CH == 1)
    Kchol = chol(Ku_s_aug);
    sol = Kchol \ (Kchol' \ RHS_s_aug);
  else
    sol = Ku_s_aug \ RHS_s_aug;
  end
  Phih_s = sol(1:ndof_tot); % Ignoring Lagrange multiplier
end

Phih(sU,1) = Phih_s(1:ndof_tot);


fprintf(' Computing the mean of Phih ((only useful for Neaumann or sanity check)).\n');

one_vec = ones(ndof_tot,1);

numerateur = Phih_s' * Mu * one_vec;
denominateur = one_vec' * Mu * one_vec;

mean_Phi_h = numerateur / denominateur;

fprintf('Average value of Phi_h on the domain = %.12e\n', mean_Phi_h);

%%%%%%%%%%%

[Estim_err, Estim_err_C, Estim_err_NC,Estim_err_data,Estim_err_without_curl] = comp_estim1(Phih,idof);
fig = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing the errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dPhi=Phih-Phiexh;


GradNorm2 = GradL2Error(Phih, idof); %% Error  H1 || Phi_exact - Phi_h || without projection on DG basis
GradNorm = sqrt(sum(GradNorm2));

Er0=sqrt(dPhi'*Mu*dPhi);
Er1=sqrt(dPhi'*Ku*dPhi);

fprintf('------------------------------------------------\n');
fprintf('P%i DG mesh_%i, ||Ph(Phiex)-Phih||_0= %7.2e\n',ordre,mi,Er0);
fprintf('P%i DG mesh_%i, ||Ph(Phiex)-Phih||_h = %7.2e\n',ordre,mi,Er1);
fprintf('P%i DG mesh_%i, sqrt( sum_T ||Grad(Ph(Phiex)-Phih)||_T^2 )= %7.2e\n',ordre,mi,GradNorm);
fprintf('P%i DG mesh_%i, Estim_err= %7.2e\n',ordre,mi,Estim_err);
fprintf('------------------------------------------------\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vizualisation of the numerical results between Phi_ex (Lagrange), Phi_h (DG) and the exact solution evaluated directly on the mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global visu
if (visu == 2 && fig > 0)
  fig+=1;
  %% Mass matrix LG
  MLG = MatMassLG(ordre);
  Phiex_LG = DGXYtoLG(Phiexh,ordre,OrdTri,MLG,idof);
  Phih_LG=DGXYtoLG(Phih,ordre,OrdTri,MLG,idof);
  if (ordre==1)
    global NumEdg
    NT=NumTri; CN=CoorNeu;
  end
  if (ordre==2)
    NT=NumTri2; CN=CoorNeu2;
  end

  tex=sprintf('Exact solution Phi, P%i, mesh%i', ordre,mi);
  tDG=sprintf('Exact P%i DG Phi, mesh%i',ordre,mi);
  tSol = sprintf('Target solution');
  figure(fig)
  subplot(1,3,1)
  colormap ("jet");
  trisurf(NT,CN(:,1),CN(:,2),Phiex_LG);
  view(2);
  shading interp
  title(tex)
  colorbar;
  %
  subplot(1,3,2)
  trisurf(NT,CN(:,1),CN(:,2),Phih_LG);
  view(2);
  shading interp
  title(tDG)
  colorbar;
  %
  subplot(1,3,3)
  trisurf(NT,CN(:,1),CN(:,2),evalSolution(CN(:,1),CN(:,2)));
  view(2);
  shading interp
  title(tSol)
  colorbar;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


