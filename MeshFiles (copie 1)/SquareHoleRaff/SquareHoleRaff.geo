// Paramètres de maillage
h = 0.1;  // Taille globale
a = 10;   // Raffinement local (facteur)

// Paramètres supplémentaires pour gérer plusieurs zones de raffinement
// Exemple : deux zones distinctes (autour du trou + autre endroit)
SetFactory("OpenCASCADE");

// Points du carré extérieur
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};

// Points du carré intérieur
Point(5) = {0.4, 0.4, 0, h/a};
Point(6) = {0.6, 0.4, 0, h/a};
Point(7) = {0.6, 0.6, 0, h/a};
Point(8) = {0.4, 0.6, 0, h/a};

// Exemple de NOUVEAUX points de raffinement (à toi de les placer où tu veux)
Point(20) = {0.15, 0.85, 0, h/a};  
Point(21) = {0.85, 0.85, 0, h/a}; 
Point(22) = {0.5, 0.8, 0, h/a};
Point(23) = {0.85, 0.6, 0, h/a};
Point(24) = {0.15, 0.6, 0, h/a};
Point(25) = {0.6, 0.6, 0, h/a};
Point(26) = {0.4, 0.6, 0, h/a};
Point(27) = {0.05, 0.95, 0, h/a};
Point(28) = {0.35, 0.95, 0, h/a};
Point(29) = {0.55, 0.95, 0, h/a};
Point(30) = {0.75, 0.95, 0, h/a};
Point(31) = {0.85, 0.15, 0, h/a};
Point(32) = {0.15, 0.15, 0, h/a};
Point(33) = {0.95, 0.95, 0, h/a};

// Lignes extérieures
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(10) = {3, 4};
Line(4) = {4, 1};

// Lignes trou intérieur
Line(11) = {5, 6};
Line(12) = {6, 7};
Line(13) = {7, 8};
Line(14) = {8, 5};

// Boucles fermées
Line Loop(1) = {1, 2, 10, 4};     // extérieur
Line Loop(2) = {11, 12, 13, 14};  // trou

// Surface avec trou
Plane Surface(1) = {1, 2};

// ====================
// Champs de raffinement
// ====================

// Zone 1 : autour du carré intérieur
Field[1] = Distance;
Field[1].PointsList = {3, 4, 5, 6, 7, 8};
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h/(3*a);
Field[2].SizeMax = h;
Field[2].DistMin = 0.05;
Field[2].DistMax = 0.15;

// Zone 2 : autour d’un point libre (Point 20)
Field[3] = Distance;
Field[3].PointsList = {20,21,22,23,24,25,26,27,28,29,30,31,32,33};
Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = h/a;
Field[4].SizeMax = h;
Field[4].DistMin = 0.1;
Field[4].DistMax = 0.5;

// Combinaison des champs
Field[100] = Min;
Field[100].FieldsList = {2,4,6};

Background Field = 100;

// ====================
// Groupes physiques
// ====================
Physical Line("Bottom") = {1};
Physical Line("Right") = {2};
Physical Line("Top") = {10};
Physical Line("Left") = {4};
Physical Line("Inner_Bottom") = {11};
Physical Line("Inner_Right") = {12};
Physical Line("Inner_Top") = {13};
Physical Line("Inner_Left") = {14};
Physical Surface("Domain") = {1};

// ====================
// Options globales
// ====================
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;

