
lc = 0.2;


Function Ball

p1 = newp; Point(p1) = {x,y,z,lc};

p2 = newp; Point(p2) = {x+r,y,z,lc};
p3 = newp; Point(p3) = {x,y+r,z,lc};
p4 = newp; Point(p4) = {x,y,z+r,lc};
p5 = newp; Point(p5) = {x-r,y,z,lc};
p6 = newp; Point(p6) = {x,y-r,z,lc};
p7 = newp; Point(p7) = {x,y,z-r,lc};

c1 = newreg; Circle(c1) = {p2,p1,p7};
c2 = newreg; Circle(c2) = {p7,p1,p5};
c3 = newreg; Circle(c3) = {p5,p1,p4}; 
c4 = newreg; Circle(c4) = {p4,p1,p2};
c5 = newreg; Circle(c5) = {p2,p1,p3}; 
c6 = newreg; Circle(c6) = {p3,p1,p5};
c7 = newreg; Circle(c7) = {p5,p1,p6}; 
c8 = newreg; Circle(c8) = {p6,p1,p2};
c9 = newreg; Circle(c9) = {p7,p1,p3}; 
c10 = newreg; Circle(c10) = {p3,p1,p4};
c11 = newreg; Circle(c11) = {p4,p1,p6}; 
c12 = newreg; Circle(c12) = {p6,p1,p7};


l1 = newl; Line Loop(l1) = {c5,c10,c4};    
s1 = news; Surface(s1) = {l1};
l2 = newl; Line Loop(l2) = {c9,-c5,c1};    
s2 = news; Surface(s2) = {l2};
l3 = newl; Line Loop(l3) = {c12,-c8,-c1}; 
s3 = news; Surface(s3) = {l3};
l4 = newl; Line Loop(l4) = {c8,-c4,c11};  
s4 = news; Surface(s4) = {l4};
l5 = newl; Line Loop(l5) = {-c10,c6,c3};  
s5 = news; Surface(s5) = {l5};
l6 = newl; Line Loop(l6) = {-c11,-c3,c7};  
s6 = news; Surface(s6) = {l6};
l7 = newl; Line Loop(l7) = {-c2,-c7,-c12};
s7 = news; Surface(s7) = {l7};
l8 = newl; Line Loop(l8) = {-c6,-c9,c2};  
s8 = news; Surface(s8) = {l8};

v1 = newv;
Return


x = 0 ; y = 0 ; z = 0 ; r = 2 ;

Call Ball;
theloops[0] = newreg;
Surface Loop(theloops[0]) = {s8,s5,s1,s2,s3,s7,s6,s4};
Printf("(center = {%g,%g,%g}, radius = %g) has number %g!", x, y, z, r, v1) ;

Volume(1) = {theloops[0]};

Physical Volume(1) = v1; 
