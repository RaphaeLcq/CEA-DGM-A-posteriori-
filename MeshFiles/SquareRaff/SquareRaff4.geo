a = 3;
h = 0.05;

// Laisser les deux premières lignes vides.
// h=0.1;
// h=0.05;
// h=0.025;
// h=0.0125;
// h=0.00625;

SetFactory("OpenCASCADE");

// Coins
Point(1) = {0, 0, 0, h/a};
Point(2) = {1, 0, 0, h/a};
Point(3) = {1, 1, 0, h/a};
Point(4) = {0, 1, 0, h/a};

// Milieux bords
Point(5) = {0, 0.5, 0, h/a};
Point(6) = {0.5, 0, 0, h/a};
Point(7) = {1, 0.5, 0, h/a};
Point(8) = {0.5, 1, 0, h/a};

// Centre
Point(9) = {0.5, 0.5, 0, h/a};

// 1/3 et 2/3
Point(10) = {1/3, 1/3, 0, h};
Point(11) = {1/3, 2/3, 0, h};
Point(12) = {2/3, 1/3, 0, h};
Point(13) = {2/3, 2/3, 0, h};

// Cercle moyen
Point(14) = {0.5 + 2/6 * Cos(0*2*Pi/9), 0.5 + 2/6 * Sin(0*2*Pi/9), 0, h};
Point(15) = {0.5 + 2/6 * Cos(1*2*Pi/9), 0.5 + 2/6 * Sin(1*2*Pi/9), 0, h};
Point(16) = {0.5 + 2/6 * Cos(2*2*Pi/9), 0.5 + 2/6 * Sin(2*2*Pi/9), 0, h};
Point(17) = {0.5 + 2/6 * Cos(3*2*Pi/9), 0.5 + 2/6 * Sin(3*2*Pi/9), 0, h};
Point(18) = {0.5 + 2/6 * Cos(4*2*Pi/9), 0.5 + 2/6 * Sin(4*2*Pi/9), 0, h};
Point(19) = {0.5 + 2/6 * Cos(5*2*Pi/9), 0.5 + 2/6 * Sin(5*2*Pi/9), 0, h};
Point(20) = {0.5 + 2/6 * Cos(6*2*Pi/9), 0.5 + 2/6 * Sin(6*2*Pi/9), 0, h};
Point(21) = {0.5 + 2/6 * Cos(7*2*Pi/9), 0.5 + 2/6 * Sin(7*2*Pi/9), 0, h};
Point(22) = {0.5 + 2/6 * Cos(8*2*Pi/9), 0.5 + 2/6 * Sin(8*2*Pi/9), 0, h};

// Cercle grand
Point(23) = {0.5 + 0.45 * Cos(0*2*Pi/12), 0.5 + 0.45 * Sin(0*2*Pi/12), 0, h};
Point(24) = {0.5 + 0.45 * Cos(1*2*Pi/12), 0.5 + 0.45 * Sin(1*2*Pi/12), 0, h};
Point(25) = {0.5 + 0.45 * Cos(2*2*Pi/12), 0.5 + 0.45 * Sin(2*2*Pi/12), 0, h};
Point(26) = {0.5 + 0.45 * Cos(3*2*Pi/12), 0.5 + 0.45 * Sin(3*2*Pi/12), 0, h};
Point(27) = {0.5 + 0.45 * Cos(4*2*Pi/12), 0.5 + 0.45 * Sin(4*2*Pi/12), 0, h};
Point(28) = {0.5 + 0.45 * Cos(5*2*Pi/12), 0.5 + 0.45 * Sin(5*2*Pi/12), 0, h};
Point(29) = {0.5 + 0.45 * Cos(6*2*Pi/12), 0.5 + 0.45 * Sin(6*2*Pi/12), 0, h};
Point(30) = {0.5 + 0.45 * Cos(7*2*Pi/12), 0.5 + 0.45 * Sin(7*2*Pi/12), 0, h};
Point(31) = {0.5 + 0.45 * Cos(8*2*Pi/12), 0.5 + 0.45 * Sin(8*2*Pi/12), 0, h};
Point(32) = {0.5 + 0.45 * Cos(9*2*Pi/12), 0.5 + 0.45 * Sin(9*2*Pi/12), 0, h};
Point(33) = {0.5 + 0.45 * Cos(10*2*Pi/12), 0.5 + 0.45 * Sin(10*2*Pi/12), 0, h};
Point(34) = {0.5 + 0.45 * Cos(11*2*Pi/12), 0.5 + 0.45 * Sin(11*2*Pi/12), 0, h};

// Cercle petit
Point(35) = {0.5 + 0.15 * Cos(0*2*Pi/8), 0.5 + 0.15 * Sin(0*2*Pi/8), 0, h};
Point(36) = {0.5 + 0.15 * Cos(1*2*Pi/8), 0.5 + 0.15 * Sin(1*2*Pi/8), 0, h};
Point(37) = {0.5 + 0.15 * Cos(2*2*Pi/8), 0.5 + 0.15 * Sin(2*2*Pi/8), 0, h};
Point(38) = {0.5 + 0.15 * Cos(3*2*Pi/8), 0.5 + 0.15 * Sin(3*2*Pi/8), 0, h};
Point(39) = {0.5 + 0.15 * Cos(4*2*Pi/8), 0.5 + 0.15 * Sin(4*2*Pi/8), 0, h};
Point(40) = {0.5 + 0.15 * Cos(5*2*Pi/8), 0.5 + 0.15 * Sin(5*2*Pi/8), 0, h};
Point(41) = {0.5 + 0.15 * Cos(6*2*Pi/8), 0.5 + 0.15 * Sin(6*2*Pi/8), 0, h};
Point(42) = {0.5 + 0.15 * Cos(7*2*Pi/8), 0.5 + 0.15 * Sin(7*2*Pi/8), 0, h};

// Contour carré
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// ----- RAFFINEMENTS -----

// Distance et seuil pour TOUS les points
Field[1] = Distance;
Field[1].PointsList = {1:42};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h/(a*2);
Field[2].SizeMax = h;
Field[2].DistMin = 0.02;
Field[2].DistMax = 0.3;

// Appliquer
Background Field = 2;

Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;

Physical Surface("Domain") = {1};

