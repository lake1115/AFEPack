MyPoint=1;
MyLine=1;
MySurface=1;
MyVolume=-1;

lc = 16;
lcc = 4;
r = 25;
R = 100;
x_center = 0;
y_center = 0;
z_center = 0;

Point(1) = {x_center, y_center, z_center, lc};
Point(2) = {x_center-R, y_center, z_center, lc};
Point(4) = {x_center, y_center-R, z_center, lc};
Point(5) = {x_center + R, y_center, z_center, lc};
Point(8) = {x_center, y_center, z_center-R, lc};
Point(11) = {x_center, y_center + R, z_center, lc};
Point(14) = {x_center, y_center, z_center + R, lc};

Point(101) = {x_center, y_center, z_center, lcc};
Point(102) = {x_center - r, y_center, z_center, lcc};
Point(104) = {x_center, y_center - r, z_center, lcc};
Point(105) = {x_center + r, y_center, z_center, lcc};
Point(108) = {x_center, y_center, z_center - r, lcc};
Point(111) = {x_center, y_center + r, z_center, lcc};
Point(114) = {x_center, y_center, z_center + r, lcc};

Physical Point(MyPoint)={1,2,4,5,8,11,14,101,102,104,105,108,111,114};

Circle (1) = {2, 1, 4} Plane{0, 0, 1};
Circle (2) = {4, 1, 5} Plane{0, 0, 1};
Circle (3) = {2, 1, 8} Plane{0, 0, 1};
Circle (4) = {4, 1, 8} Plane{0, 0, 1};
Circle (6) = {2, 1, 11} Plane{0, 0, 1};
Circle (7) = {8, 1, 11} Plane{0, 0, 1};
Circle (9) = {2, 1, 14} Plane{0, 0, 1};
Circle (10) = {11, 1, 14} Plane{0, 0, 1};
Circle (13) = {14, 1, 4} Plane{0, 0, 1};
Circle (15) = {8, 1, 5} Plane{0, 0, 1};
Circle (18) = {11, 1, 5} Plane{0, 0, 1};
Circle (21) = {14, 1, 5} Plane{0, 0, 1};

Circle (101) = {102, 101, 104} Plane{0, 0, 1};
Circle (102) = {104, 101, 105} Plane{0, 0, 1};
Circle (103) = {102, 101, 108} Plane{0, 0, 1};
Circle (104) = {104, 101, 108} Plane{0, 0, 1};
Circle (106) = {102, 101, 111} Plane{0, 0, 1};
Circle (107) = {108, 101, 111} Plane{0, 0, 1};
Circle (109) = {102, 101, 114} Plane{0, 0, 1};
Circle (110) = {111, 101, 114} Plane{0, 0, 1};
Circle (113) = {114, 101, 104} Plane{0, 0, 1};
Circle (115) = {108, 101, 105} Plane{0, 0, 1};
Circle (118) = {111, 101, 105} Plane{0, 0, 1};
Circle (121) = {114, 101, 105} Plane{0, 0, 1};

Physical Line(MyLine) = {1,2,3,4,6,7,9,10,13,15,18,21,101,102,103,104,106,107,109,110,113,115,118,121};


Line Loop (1000005) = {1, 4, -3};
Ruled Surface (5) = {1000005};
Line Loop (1000008) = {3, 7, -6};
Ruled Surface (8) = {1000008};
Line Loop (1000011) = {6, 10, -9};
Ruled Surface (11) = {1000011};
Line Loop (1000014) = {9, 13, -1};
Ruled Surface (14) = {1000014};
Line Loop (1000017) = {-15, -4, 2};
Ruled Surface (17) = {1000017};
Line Loop (1000020) = {-18, -7, 15};
Ruled Surface (20) = {1000020};
Line Loop (1000023) = {-21, -10, 18};
Ruled Surface (23) = {1000023};
Line Loop (1000026) = {-2, -13, 21};
Ruled Surface (26) = {1000026};

Line Loop (11000005) = {101, 104, -103};
Ruled Surface (105) = {11000005};
Line Loop (11000008) = {103, 107, -106};
Ruled Surface (108) = {11000008};
Line Loop (11000011) = {106, 110, -109};
Ruled Surface (111) = {11000011};
Line Loop (11000014) = {109, 113, -101};
Ruled Surface (114) = {11000014};
Line Loop (11000017) = {-115, -104, 102};
Ruled Surface (117) = {11000017};
Line Loop (11000020) = {-118, -107, 115};
Ruled Surface (120) = {11000020};
Line Loop (11000023) = {-121, -110, 118};
Ruled Surface (123) = {11000023};
Line Loop (11000026) = {-102, -113, 121};
Ruled Surface (126) = {11000026};

 

Surface Loop(1) = { 5, 8, 11, 14, 17, 20, 23, 26 };
Surface Loop(2) = { 105, 108, 111, 114, 117, 120, 123, 126 };
Volume(1) = { 1, -2};

Physical Surface(MySurface) = {23, 11, 8, 20, 17, 26, 5, 14, 114, 111, 108, 117, 105, 126, 123, 120};
Physical Volume(MyVolume) = {1};
