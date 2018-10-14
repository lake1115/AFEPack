
lc = 1e-1;


Point(1) = {-1, -1, 0, lc};
Point(2) = {1, -1, 0, lc} ;
Point(3) = {1, 1, 0, lc} ;
Point(4) = {-1, 1, 0, lc} ;

Point(5) = {-0.5, -0.5, 0, lc};
Point(6) = {0.5, -0.5, 0, lc} ;
Point(7) = {0.5, 0.5, 0, lc} ;
Point(8) = {-0.5, 0.5, 0, lc} ;


Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

Line(5) = {5,6} ;
Line(6) = {6,7} ;
Line(7) = {7,8} ;
Line(8) = {8,5} ;

Curve Loop(1) = {1,2,3,4} ;

Curve Loop(2) = {5,6,7,8} ;

Plane Surface(1) = {1,2} ;

