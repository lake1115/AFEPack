MyPoint=1;
MyLine=1;
MySurface=1;
MyVolume=-1;
lc=0.2;

Point(1) = {-1, -1, -1, lc};
Point(2) = { 1, -1, -1, lc};
Point(3) = { 1,  1, -1, lc};
Point(4) = {-1,  1, -1, lc};
Point(5) = {-1, -1,  1, lc};
Point(6) = { 1, -1,  1, lc};
Point(7) = { 1,  1,  1, lc};
Point(8) = {-1,  1,  1, lc};

Physical Point(MyPoint) = {1,2,3,4,5,6,7,8};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3)={3,4};
Line(4)={4,1};
Line(5) = {1,5};
Line(6) = {5,6};
Line(7)= {6,7};
Line(8)= {7,8};
Line(9) = {8,5};
Line(10) = {2,6};
Line(11)= {3,7};
Line(12)= {4,8};

Physical Line(MyLine)= {1,2,3,4,5,6,7,8,9,10,11,12};

Line Loop(13) = {1,2,3,4};
Plane Surface(14) = {13};
Line Loop(15) = {1,10,-6,-5};
Plane Surface(16) = {15};
Line Loop(17) = {2,11,-7,-10};
Plane Surface(18) = {17};
Line Loop(19) = {3,12,-8,-11};
Plane Surface(20) = {19};
Line Loop(21) = {4,5,-9,-12};
Plane Surface(22) = {21};
Line Loop(23) = {6,7,8,9};
Plane Surface(24) = {23};

Physical Surface(MySurface) = {20,14,16,18,24,22};

Surface Loop(25) = {20,14,16,18,24,22};
Volume(26) = {25};

Physical Volume(MyVolume) = {26};
