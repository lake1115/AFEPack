T_Point_bmark = 2;
D_Point_bmark = 1;
T_Line_bmark = 2;
D_Line_bmark = 1;
TBC_bmark = 2;
Dirichlet_bmark = 1;

lc = 0.2;
r = 0.5;
R = 1;
h = 1;

//Point
//cylinder
Point(1) = {0,0,0,lc};

Point(2) = {r,0,0,lc};
Point(3) = {0,r,0,lc};
Point(4) = {-r,0,0,lc};
Point(5) = {0,-r,0,lc};

Point(6) = {0,0,-h,lc};

Point(7) = {r,0,-h,lc};
Point(8) = {0,r,-h,lc};
Point(9) = {-r,0,-h,lc};
Point(10) = {0,-r,-h,lc};
//half-ball
Point(11) = {R,0,0,lc};
Point(12) = {0,R,0,lc};
Point(13) = {-R,0,0,lc};
Point(14) = {0,-R,0,lc};
Point(15) = {0,0,R,lc};

Physical Point(T_Point_bmark) = {11,12,13,14,15};
Physical Point(D_Point_bmark) = {2,3,4,5,7,8,9,10};
//Line
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Circle(5) = {7,6,8};
Circle(6) = {8,6,9};
Circle(7) = {9,6,10};
Circle(8) = {10,6,7};

Line(9) = {2,7};
Line(10) = {3,8};
Line(11) = {4,9};
Line(12) = {5,10};

Circle(13) = {11,1,12};
Circle(14) = {12,1,13};
Circle(15) = {13,1,14};
Circle(16) = {14,1,11};
Circle(17) = {11,1,15};
Circle(18) = {12,1,15};
Circle(19) = {13,1,15};
Circle(20) = {14,1,15};

Line(21) = {11,2};
Line(22) = {12,3};
Line(23) = {13,4};
Line(24) = {14,5};

Physical Line(T_Line_bmark) = {13,14,15,16,17,18,19,20};
Physical Line(D_Line_bmark) = {1,2,3,4,5,6,7,8,9,10,11,12,21,22,23,24};
//Surface
Line Loop(25) = {9,5,-10,-1};
Surface(26) = {25};
Line Loop(27) = {10,6,-11,-2};
Surface(28) = {27};
Line Loop(29) = {11,7,-12,-3};
Surface(30) = {29}; 
Line Loop(31) = {12,8,-9,-4};
Surface(32) = {31};
Line Loop(33) = {5,6,7,8};
Surface(34) = {33};

Line Loop(35) = {21,1,-22,-13};
Surface(36) = {35};
Line Loop(37) = {22,2,-23,-14};
Surface(38) = {37};
Line Loop(39) = {23,3,-24,-15};
Surface(40) = {39};
Line Loop(41) = {24,4,-21,-16};
Surface(42) = {41};

Line Loop(43) = {13,18,-17};
Surface(44) = {43};
Line Loop(45) = {14,19,-18};
Surface(46) = {45};
Line Loop(47) = {15,20,-19};
Surface(48) = {47};
Line Loop(49) = {16,17,-20};
Surface(50) = {49};

Physical Surface(TBC_bmark) = {44,46,48,50};
Physical Surface(Dirichlet_bmark) = {34,26,28,30,32,36,38,40,42}; 

//Volume
Surface Loop(51) = {34,26,28,30,32,36,38,40,42,44,46,48,50};
Volume(52) = {51};
Physical Volume(1) = {52};

