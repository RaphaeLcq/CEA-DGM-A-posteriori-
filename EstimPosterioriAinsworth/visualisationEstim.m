

function visualisationEstim(visu, Estim_err_NC2, curl_rhok_norm2, norm_rho_K2, contribution_rhok2, norm2_Data, Estim_err_C2_without_curl, Estim_err_C2, Estim_err_tot2, BrokenNorm2Real)

global NumTri CoorNeu

if visu == 1 %% Vizualisation


  %%%% Visualisation des différentes contribution de l'erreur %%%%%%%%%%%%%%%%%%%
  fig = 100;
  figure(fig);

  subplot(2,3,1)
  patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', sqrt(norm_rho_K2), 'FaceColor', 'flat', 'EdgeColor', 'k');
  axis equal;
  view(2);
  colormap('jet');
  colorbar;
  title('Carte locale de norme L2 de rho_K');


  subplot(2,3,2)
  patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', sqrt(curl_rhok_norm2), 'FaceColor', 'flat', 'EdgeColor', 'k');
  axis equal;
  view(2);
  colormap('jet');
  colorbar;
  title('Carte locale de norme L2 de C_t * |K| * curl(rho_K)');


  subplot(2,3,3)
  patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', sqrt(contribution_rhok2), 'FaceColor', 'flat', 'EdgeColor', 'k');
  axis equal;
  view(2);
  colormap('jet');
  colorbar;
  title('Carte locale de la contribution totale associée à rho_K');

  subplot(2,3,4)
  patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', sqrt(norm2_Data), 'FaceColor', 'flat', 'EdgeColor', 'k');
  axis equal;
  view(2);
  colormap('jet');
  colorbar;
  title('Carte locale de l''oscillation des données');

  subplot(2,3,5)
  patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', sqrt(Estim_err_C2_without_curl), 'FaceColor', 'flat', 'EdgeColor', 'k');
  axis equal;
  view(2);
  colormap('jet');
  colorbar;
  title('Carte locale de l''estimateur sans l''amélioration avec le rotationnel');

  subplot(2,3,6)
  patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', sqrt(Estim_err_C2), 'FaceColor', 'flat', 'EdgeColor', 'k');
  axis equal;
  view(2);
  colormap('jet');
  colorbar;
  title('Carte locale de l''estimateur avec l''amélioration avec le rotationnel');

  fig=fig-1;
  figure(fig);
  subplot(1,2,1)
  patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', sqrt(Estim_err_NC2), 'FaceColor', 'flat', 'EdgeColor', 'k');
  axis equal;
  view(2);
  colormap('jet');
  colorbar;
  title('Carte locale de l’estimateur d''erreur non conforme');


  subplot(1,2,2)
  patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', sqrt(Estim_err_C2), 'FaceColor', 'flat', 'EdgeColor', 'k');
  axis equal;
  view(2);
  colormap('jet');
  colorbar;
  title('Carte locale de l''estimation d''erreur conforme');


  %%%% Visualisation de l'erreur réelle du gradient brisé et de l'estimation de l'erreur
  min_val = min([sqrt(Estim_err_tot2(:)); sqrt(BrokenNorm2Real(:))]);
  max_val = max([sqrt(Estim_err_tot2(:)); sqrt(BrokenNorm2Real(:))]);


  fig=fig-1;
  figure(fig);
  subplot(1,2,1)
  patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', sqrt(Estim_err_tot2), 'FaceColor', 'flat', 'EdgeColor', 'k');
  axis equal;
  view(2);
  colormap('jet');
##  caxis([min_val max_val]); % Appliquer la même échelle de couleurs
  colorbar;
  title('Carte locale de l''estimation d''erreur totale');

  subplot(1,2,2)
  patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', sqrt(BrokenNorm2Real), 'FaceColor', 'flat', 'EdgeColor', 'k');
  axis equal;
  view(2);
  colormap('jet');
##  caxis([min_val max_val]); % Appliquer la même échelle de couleurs
  colorbar;
  title('Err grad avec la sol exacte');

  fig = fig -1;
  figure(fig);
    % Calcul des moyennes pour chaque jeu de données
  mean_err = mean(sqrt(Estim_err_tot2(:)));
  mean_broken = mean(sqrt(BrokenNorm2Real(:)));

  % Seuillage des données (mettre à 0 si <= moyenne)
  thresholded_err = sqrt(Estim_err_tot2);
  thresholded_err(thresholded_err <= mean_err) = 0;

  thresholded_broken = sqrt(BrokenNorm2Real);
  thresholded_broken(thresholded_broken <= mean_broken) = 0;

  % Trouver les valeurs max pour l'échelle de couleur commune
  max_val = max([thresholded_err(:); thresholded_broken(:)]);

  subplot(1,2,1)
  patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', thresholded_err, 'FaceColor', 'flat', 'EdgeColor', 'k');
  axis equal;
  view(2);
  colormap('jet');
  caxis([0 max_val]); % Échelle commune avec 0 comme minimum
  colorbar;
  title('Erreurs > moyenne (estimation totale)');

  subplot(1,2,2)
  patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', thresholded_broken, 'FaceColor', 'flat', 'EdgeColor', 'k');
  axis equal;
  view(2);
  colormap('jet');
  caxis([0 max_val]); % Échelle commune avec 0 comme minimum
  colorbar;
  title('Erreurs > moyenne (gradient exact)');
  global Lshape
  if (Lshape == 1) %%%% Pour visualiser sans l'erreur dans le coin
    thresholded_errLshape = sqrt(Estim_err_tot2);
    thresholded_errLshape(thresholded_errLshape >= 2*mean_err) = 0;

    thresholded_brokenLshape = sqrt(BrokenNorm2Real);
    thresholded_brokenLshape(thresholded_brokenLshape >= 2*mean_broken) = 0;

    fig=fig-1;
    figure(fig);
    subplot(1,2,1)
    patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', thresholded_errLshape, 'FaceColor', 'flat', 'EdgeColor', 'k');
    axis equal;
    view(2);
    colormap('jet');
  ##  caxis([min_val max_val]); % Appliquer la même échelle de couleurs
    colorbar;
    title('Carte locale de l''estimation d''erreur totale sans le coin');

    subplot(1,2,2)
    patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', thresholded_brokenLshape, 'FaceColor', 'flat', 'EdgeColor', 'k');
    axis equal;
    view(2);
    colormap('jet');
  ##  caxis([min_val max_val]); % Appliquer la même échelle de couleurs
    colorbar;
    title('Err grad avec la sol exacte sans le coin');
  endif

  fig=fig-1;
  figure(fig);
  patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceColor', 'white', 'EdgeColor', 'k');
  axis equal;
  view(2);
  title('Maillage triangulaire');




endif

endfunction

