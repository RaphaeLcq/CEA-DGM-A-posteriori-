// Paramètres de maillage
h = 0.00625;  // Taille de maille globale (fixe)


//h=0.1;
//h=0.05;
//h=0.025;
//h=0.0125;
//h=0.00625; 
// Paramètres de maillage



// Usine de géométrie
SetFactory("OpenCASCADE");

// Points
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};

// Lignes (numérotation inchangée)
Line(10) = {1, 2};
Line(11) = {2, 3};
Line(12) = {3, 4};
Line(13) = {4, 1};

// Boucle et surface
Line Loop(1) = {10, 11, 12, 13};
Plane Surface(1) = {1};

// Définir la taille min et max pour le maillage global
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;

// Groupes physiques
Physical Line("Bottom") = {10};
Physical Line("Right")  = {11};
Physical Line("Top")    = {12};
Physical Line("Left")   = {13};
Physical Surface("Domain") = {1};
