// Paramètres de maillage
h = 0.00625;  // Taille de maille globale (fixe)
a = 1;     // SEUL ce paramètre contrôle le raffinement local

//h=0.1;
//h=0.05;
//h=0.025;
//h=0.0125;
//h=0.00625; 

SetFactory("OpenCASCADE");

// Points du carré extérieur (taille fixe h)
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};

// Points du carré intérieur (taille h/a - SEULS CES POINTS SONT RAFFINÉS)
Point(5) = {0.4, 0.4, 0, h/a};
Point(6) = {0.6, 0.4, 0, h/a};
Point(7) = {0.6, 0.6, 0, h/a};
Point(8) = {0.4, 0.6, 0, h/a};

// Lignes du contour extérieur
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(10) = {3, 4};
Line(4) = {4, 1};

// Lignes du trou
Line(11) = {5, 6};
Line(12) = {6, 7};
Line(13) = {7, 8};
Line(14) = {8, 5};

// Boucles fermées
Line Loop(1) = {1, 2, 10, 4};      // Contour extérieur
Line Loop(2) = {11, 12, 13, 14};      // Trou

// Surface avec trou
Plane Surface(1) = {1, 2};

// Champ de raffinement LOCAL uniquement autour des points 5-8
Field[1] = Distance;
Field[1].PointsList = {5, 6, 7, 8};  // Raffinement SEULEMENT autour de ces points

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h/a;        // Taille près du trou (contrôlée par 'a')
Field[2].SizeMax = h;          // Taille loin du trou (fixe)
Field[2].DistMin = 0.02;       // Rayon de raffinement maximal
Field[2].DistMax = 0.15;       // Distance de transition

Background Field = 2;

// Groupes physiques
Physical Line("Bottom") = {1};
Physical Line("Right") = {2};
Physical Line("Top") = {10};
Physical Line("Left") = {4};
Physical Line("Inner_Bottom") = {11};
Physical Line("Inner_Right") = {12};
Physical Line("Inner_Top") = {13};
Physical Line("Inner_Left") = {14};
Physical Surface("Domain") = {1};

// Options pour assurer le contrôle local
Mesh.CharacteristicLengthFromPoints = 0;  // Ignore les tailles aux points
Mesh.CharacteristicLengthExtendFromBoundary = 0;
