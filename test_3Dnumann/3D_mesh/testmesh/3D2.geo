Point_D=1;
Point_N=2;
Line_D = 1;
Line_N = 2;
Surface_D=1;
Surface_N = 2;
Volume_bmark=-1;


lc = 0.5;


Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {0, 1, 0, lc};
Point(4) = {0, 0, 1, lc};

Physical Point(Point_D) = {1,2,3,4};

Line(1) = {1,2};
Line(2) = {1,3};
Line(3) = {1,4};
Line(4) = {2,3};
Line(5) = {2,4};
Line(6) = {3,4};

Physical Line(Line_D)= {1,2,3,4,5,6};

Line Loop(1) = {3,-5,-1} ;
Plane Surface(1) = {1} ;
Line Loop(2) = {2,6,-3};
Plane Surface(2) = {2};
Line Loop(3) = {1,4,-2};
Plane Surface(3) = {3};
Line Loop(4) = {4,6,-5};
Plane Surface(4) = {4};


Physical Surface(Surface_D) = {1,2,3,4};

Surface Loop(1) = {1,2,3,4};

Volume(1) = {1};

Physical Volume(Volume_bmark) = {1};
