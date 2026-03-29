h=0.1;
a=600;

//h=0.1;
//h=0.05;
//h=0.025;
//h=0.0125;
//h=0.00625;

SetFactory("OpenCASCADE");

// ====================
// Points
// ====================
Point(1)={0, 0, 0, h};
Point(2)={1/2, 0, 0, h};
Point(3)={1/2, 1/2, 0, h};   // coin rentrant
Point(4)={1, 1/2, 0, h};
Point(5)={1, 1, 0, h};
Point(6)={0, 1, 0, h};

// ====================
// Lignes et surface
// ====================
Line(10) ={1,2};
Line(11) ={2,3};
Line(12) ={3,4};
Line(13) ={4,5};
Line(14) ={5,6};
Line(15) ={6,1};

Line Loop(16)={10,11,12,13,14,15};
Plane Surface(17)= {16}; 

// ====================
// Raffinement local au coin (0.5,0.5)
// ====================
Field[1] = Distance;
Field[1].PointsList = {3};        // raffinement centré sur le coin rentrant

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h/a;           // taille très fine
Field[2].SizeMax = h;             // taille globale ailleurs
Field[2].DistMin = 0.01;          // rayon de raffinement
Field[2].DistMax = 0.55;          // transition douce

Background Field = 2;

// ====================
// Groupes physiques
// ====================
Physical Line("BD")={10};
Physical Line("BR")={11};
Physical Line("LD")={12};
Physical Line("LR")={13};
Physical Line("BU")={14};
Physical Line("BL")={15};
Physical Surface("Domain")={17};

// ====================
// Options globales
// ====================
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;

