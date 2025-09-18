h = 0.05;  // Modifiez cette valeur


//h=0.1;
//h=0.05;
//h=0.025;
//h=0.0125;
//h=0.00625;


SetFactory("OpenCASCADE");

Point(1)={0, 0, 0, h};
Point(2)={1/2, 0, 0, h};
Point(3)={1/2, 1/2, 0, h};
Point(4)={1, 1/2, 0, h};
Point(5)={1, 1, 0, h};
Point(6)={0, 1, 0, h};

Line(10) ={1,2};
Line(11) ={2,3};
Line(12) ={3,4};
Line(13) ={4,5};
Line(14) ={5,6};
Line(15) ={6,1};

Line Loop(16)={10,11,12,13,14,15};
Plane Surface(17)= {16};

Physical Line("BD")={10};
Physical Line("BR")={11};
Physical Line("LD")={12};
Physical Line("LR")={13};
Physical Line("BU")={14};
Physical Surface("Domain") = {17};

// Forcer les paramètres de maillage
Mesh.CharacteristicLengthMax = h;
Mesh.CharacteristicLengthMin = h/10;
