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
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(10) = {3, 4};
Line(3) = {4, 1};

// Boucle et surface
Line Loop(1) = {1, 2, 10, 3};
Plane Surface(1) = {1};

// Définir la taille min et max pour le maillage global
Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;

// Groupes physiques
Physical Line("Bottom") = {1};
Physical Line("Right")  = {2};
Physical Line("Top")    = {10};
Physical Line("Left")   = {3};
Physical Surface("Domain") = {1};
