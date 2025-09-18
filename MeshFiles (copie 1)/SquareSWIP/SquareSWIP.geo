h = 0.0125;
// Laisser les deux première lignes vides.
//h=0.1;
//h=0.05;
//h=0.025;
//h=0.0125;
//h=0.00625;
//+
SetFactory("OpenCASCADE");
Mesh.CharacteristicLengthMax = h;
Mesh.CharacteristicLengthMin = h;
// Définition des points aux coins du carré
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};

// Lignes reliant les points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Boucle fermée
Line Loop(5) = {1, 2, 3, 4};

// Surface plane
Plane Surface(6) = {5};

// Groupes physiques pour les côtés
Physical Line("Bottom") = {1};
Physical Line("Right")  = {2};
Physical Line("Top")    = {3};
Physical Line("Left")   = {4};
Physical Surface("Domain") = {6};

