
lc = 2;


Point(1) = {-2, -2, -2, lc};
Point(2) = { 2, -2, -2, lc};
Point(3) = { 2,  2, -2, lc};
Point(4) = {-2,  2, -2, lc};
Point(5) = {-2, -2,  2, lc};
Point(6) = { 2, -2,  2, lc};
Point(7) = { 2,  2,  2, lc};
Point(8) = {-2,  2,  2, lc};

Point(9)  = {-1, -1, -1, lc};
Point(10) = { 1, -1, -1, lc};
Point(11) = { 1,  1, -1, lc};
Point(12) = {-1,  1, -1, lc};
Point(13) = {-1, -1,  1, lc};
Point(14) = { 1, -1,  1, lc};
Point(15) = { 1,  1,  1, lc};
Point(16) = {-1,  1,  1, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {1,5};
Line(6) = {5,6};
Line(7) = {6,7};
Line(8) = {7,8};
Line(9) = {8,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

Line(13) = {9,10};
Line(14) = {10,11};
Line(15) = {11,12};
Line(16) = {12,9};
Line(17) = {9,13};
Line(18) = {13,14};
Line(19) = {14,15};
Line(20) = {15,16};
Line(21) = {16,13};
Line(22) = {10,14};
Line(23) = {11,15};
Line(24) = {12,16};


Line Loop(1) = {1,2,3,4} ;
Plane Surface(1) = {1} ;
Line Loop(2) = {1,10,-6,-5};
Plane Surface(2) = {2};
Line Loop(3) = {2,11,-7,-10};
Plane Surface(3) = {3};
Line Loop(4) = {3,12,-8,-11};
Plane Surface(4) = {4};
Line Loop(5) = {4,5,-9,-12};
Plane Surface(5) = {5};
Line Loop(6) = {6,7,8,9};
Plane Surface(6) = {6};

Line Loop(7) = {13,14,15,16} ;
Plane Surface(7) = {7} ;
Line Loop(8) = {13,22,-18,-17};
Plane Surface(8) = {8};
Line Loop(9) = {14,23,-19,-22};
Plane Surface(9) = {9};
Line Loop(10) = {15,24,-20,-23};
Plane Surface(10) = {10};
Line Loop(11) = {16,17,-21,-24};
Plane Surface(11) = {11};
Line Loop(12) = {18,19,20,21};
Plane Surface(12) = {12};

Surface Loop(1) = {4,1,2,3,6,5};
Surface Loop(2) = {10,7,8,9,12,11};

Volume(1) = {1,2};

