
lc = 0.5;


Macro Cube

p1 = newp; Point(p1) = {x-r, y-r, z-r, lc};
p2 = newp; Point(p2) = {x+r, y-r, z-r, lc};
p3 = newp; Point(p3) = {x+r, y+r, z-r, lc};
p4 = newp; Point(p4) = {x-r, y+r, z-r, lc};
p5 = newp; Point(p5) = {x-r, y-r, z+r, lc};
p6 = newp; Point(p6) = {x+r, y-r, z+r, lc};
p7 = newp; Point(p7) = {x+r, y+r, z+r, lc};
p8 = newp; Point(p8) = {x-r, y+r, z+r, lc};

l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p3,p4};
l4 = newl; Line(l4) = {p4,p1};
l5 = newl; Line(l5) = {p1,p5};
l6 = newl; Line(l6) = {p5,p6};
l7 = newl; Line(l7) = {p6,p7};
l8 = newl; Line(l8) = {p7,p8};
l9 = newl; Line(l9) = {p8,p5};
l10 = newl; Line(l10) = {p2,p6};
l11 = newl; Line(l11) = {p3,p7};
l12 = newl; Line(l12) = {p4,p8};


c1 = newreg; Line Loop(c1) = {l1,l2,l3,l4};
s1 = news; Plane Surface(s1) = {c1};
c2 = newreg; Line Loop(c2) = {l1,l10,-l6,-l5};
s2 = news; Plane Surface(s2) = {c2};
c3 = newreg; Line Loop(c3) = {l2,l11,-l7,-l10};
s3 = news; Plane Surface(s3) = {c3};
c4 = newreg; Line Loop(c4) = {l3,l12,-l8,-l11};
s4 = news; Plane Surface(s4) = {c4};
c5 = newreg; Line Loop(c5) = {l4,l5,-l9,-l12};
s5 = news; Plane Surface(s5) = {c5};
c6 = newreg; Line Loop(c6) = {l6,l7,l8,l9};
s6 = news; Plane Surface(s6) = {c6};

theloops[t] = newreg;
Surface Loop(theloops[t]) = {s4,s1,s2,s3,s6,s5};

v1 = newv;
//Volume(v1) = {theloops[t]};

Return

x = 0 ; y = 0 ; z = 0 ; r = 1 ;

For t In {0:1}
	
	Call Cube;
	Printf("Hole %g (center = {%g,%g,%g}, radius = %g) has number %g!",
	 t, x, y, z, r, v1) ;

	r -= 0.1;
	Physical Volume(t) = v1;
EndFor

Volume(100) = {theloops[]};
