
Point (1) = {0., 0., 0., 1.};
Point (2) = {1., 0., 0., 1.};
Point (3) = {1., 1., 0., 1.};
Point (4) = {0., 1., 0., 1.};

Line (1) = {1, 2};
Line (2) = {2, 3};
Line (3) = {3, 4};
Line (4) = {4, 1};

Curve Loop (1) = {1, 2, 3, 4};

Surface (1) = {1};

Transfinite Curve{1} = 1000;
Transfinite Curve{3} = 1000;
Transfinite Curve{2} = 1;
Transfinite Curve{4} = 1;

Transfinite Surface{1} = {1, 2, 3, 4};

Recombine Surface{1};

Physical Curve("empty")    = {1, 3};
Physical Curve("right")  = {2};
Physical Curve("left")   = {4};
Physical Surface("internal") = {1};
