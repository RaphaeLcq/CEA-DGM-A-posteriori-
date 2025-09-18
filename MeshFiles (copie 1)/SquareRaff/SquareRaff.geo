a=2;
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
Point(1) = {0, 0, 0, h/a};
Point(2) = {1, 0, 0, h/a};
Point(3) = {1, 1, 0, h/a};
Point(4) = {0, 1, 0, h/a};

// Points milieux des bords (maillage fin)
Point(5) = {0, 0.5, 0, h/a};
Point(6) = {0.5, 0, 0, h/a};
Point(7) = {1, 0.5, 0, h/a};
Point(8) = {0.5, 1, 0, h/a};

// Centre (maillage fin)
Point(9) = {0.5, 0.5, 0, h/a};

// Lignes du contour uniquement
Line(1) = {1, 6};
Line(2) = {6, 2};
Line(3) = {2, 7};
Line(4) = {7, 3};
Line(5) = {3, 8};
Line(6) = {8, 4};
Line(7) = {4, 5};
Line(8) = {5, 1};

// Définition du contour externe - UNE SEULE SURFACE
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

// Contraintes de taille pour assurer la transition
Field[1] = Distance;
Field[1].PointsList = {1, 2, 3, 4, 5, 6, 7, 8, 9};

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
