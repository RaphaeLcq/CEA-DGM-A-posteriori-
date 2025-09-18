a=1;
h=0.00625;

// Laisser les deux première lignes vides.
//h=0.1;
//h=0.05;
//h=0.025;
//h=0.0125;
//h=0.00625;
//+
SetFactory("OpenCASCADE");
//+
Point(1)={0, 0, 0, h};
//+
Point(2)={1/2, 0, 0, h};
//+
Point(3)={1/2, 1/2, 0, h/a};
//+
Point(4)={1, 1/2, 0, h};
//+
Point(5)={1, 1, 0, h};
//+
Point(6)={0, 1, 0, h};
//+
Line(1) ={1,2};
//+
Line(2) ={2,3};
//+
Line(3) ={3,4};
//+
Line(4) ={4,5};
//+
Line(5) ={5,6};
//+
Line(6) ={6,1};
//+ 
Line Loop(7)={1,2,3,4,5,6};
//+ 
Plane Surface(8)= {7}; 
//+
Physical Line("BD")={1};
//+
Physical Line("BR")={2};
//+
Physical Line("LD")={3};
//+
Physical Line("LR")={4};
//+
Physical Line("BU")={5};
