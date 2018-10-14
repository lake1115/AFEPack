
lc = 1e-1;


Point(1) = {0, 0, 0, lc};

Point(2) = {2, 0, 0, lc} ;
Point(3) = {0, 2, 0, lc} ;
Point(4) = {-2, 0, 0, lc} ;
Point(5) = {0, -2, 0, lc};

Point(6) = {1, 0, 0, lc} ;
Point(7) = {0, 1, 0, lc} ;
Point(8) = {-1, 0, 0, lc} ;
Point(9) = {0, -1, 0, lc} ;


Circle(1) = {2,1,3} ;
Circle(2) = {3,1,4} ;
Circle(3) = {4,1,5} ;
Circle(4) = {5,1,2} ;

Circle(5) = {6,1,7} ;
Circle(6) = {7,1,8} ;
Circle(7) = {8,1,9} ;
Circle(8) = {9,1,6} ;

Curve Loop(1) = {1,2,3,4} ;

Curve Loop(2) = {5,6,7,8} ;

Plane Surface(1) = {1,2} ;

