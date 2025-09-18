a=20;
h=0.1;

// Laisser les deux première lignes vides.
//h=0.1;
//h=0.05;
//h=0.025;
//h=0.0125;
//h=0.00625;
//+
SetFactory("OpenCASCADE");

// Points aux coins (maillage raffiné)
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};

// Points segment (maillage fin)
Point(5) = {0.6, 0.5, 0, h/a};
Point(6) = {0.75, 0.5, 0, h/a};
Point(7) = {0.85, 0.5, 0, h/a};
Point(8) = {1, 0.5, 0, h/a};

// Centre (maillage fin)
Point(9) = {0.5, 0.5, 0, h/a};

// Lignes du contour uniquement
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Définition du contour externe - UNE SEULE SURFACE
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Contraintes de taille pour assurer la transition
Field[1] = Distance;
Field[1].PointsList = {5, 6, 7, 8, 9};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h/a;
Field[2].SizeMax = h;
Field[2].DistMin = 0.1;
Field[2].DistMax = 0.4;

Background Field = 2;

// Désactive la taille par défaut
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;
