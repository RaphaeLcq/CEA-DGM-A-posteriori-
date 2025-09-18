h = 0.1;
h_fine = h/100;  // Taille fine pour le raffinement local
// Laisser les deux première lignes vides.
//h=0.1;
//h=0.05;
//h=0.025;
//h=0.0125;
//h=0.00625;
//+
SetFactory("OpenCASCADE");

// Définition des points aux coins du carré
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};

// Point central pour le raffinement
Point(5) = {0.5, 0.5, 0, h_fine};

// Lignes du carré unitaire
Line(1) = {1, 2};  // Bas
Line(2) = {2, 3};  // Droite
Line(3) = {3, 4};  // Haut
Line(4) = {4, 1};  // Gauche

// Boucle fermée du domaine complet
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// IMPORTANT: Ne pas utiliser Transfinite ni Recombine pour avoir des triangles !
// Laisser gmsh créer un maillage triangulaire libre

// Raffinement local autour du point central
Field[1] = Distance;
Field[1].PointsList = {5};  // Distance au point (0.5, 0.5)

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h_fine;   // Taille fine près du centre
Field[2].SizeMax = h;        // Taille normale ailleurs
Field[2].DistMin = 0.05;     // Rayon de la zone très fine (réduit)
Field[2].DistMax = 0.2;      // Rayon de la zone de transition

// Appliquer le champ de raffinement
Background Field = 2;

// Configuration pour un maillage triangulaire adaptatif
Mesh.CharacteristicLengthFromPoints = 1;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 1;

// Forcer l'algorithme de maillage triangulaire

Mesh.ElementOrder = 1;

// S'assurer qu'on a bien des triangles
Mesh.RecombineAll = 0;  // Pas de recombinaison en quadrilatères


// Champ basé sur plusieurs points pour un meilleur contrôle
Field[3] = Distance;
Field[3].PointsList = {5, 6, 7, 8, 9};

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = h_fine;
Field[4].SizeMax = h;
Field[4].DistMin = 0.03;
Field[4].DistMax = 0.15;

// Utiliser le champ avec plus de points
Background Field = 4;

// Groupes physiques pour les côtés
Physical Line("Bottom") = {1};
Physical Line("Right")  = {2};
Physical Line("Top")    = {3};
Physical Line("Left")   = {4};
Physical Surface("Domain") = {6};
