Point_bmark=1;
Line_bmark=1;
Surface_bmark=1;
Volume_bmark=-1;


lc = 0.25;


Point(1) = {-1, -1, -1, lc};
Point(2) = { 1, -1, -1, lc};
Point(3) = { 1,  1, -1, lc};
Point(4) = {-1,  1, -1, lc};
Point(5) = {-1, -1,  1, lc};
Point(6) = { 1, -1,  1, lc};
Point(7) = { 1,  1,  1, lc};
Point(8) = {-1,  1,  1, lc};
Physical Point(Point_bmark) = {1,2,3,4,5,6,7,8};

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
Physical Line(Line_bmark)= {1,2,3,4,5,6,7,8,9,10,11,12};



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

Physical Surface(Surface_bmark) = {4,1,2,3,6,5};

Surface Loop(1) = {4,1,2,3,6,5};

Volume(1) = {1};

Physical Volume(Volume_bmark) = {1};
