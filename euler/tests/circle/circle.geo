SetFactory("OpenCASCADE");

Circle (1) = {0, 0, 0, 1};
Circle (2) = {0, 0, 0, 30};

Curve Loop (1) = {1};
Curve Loop (2) = {2};

Plane Surface (1) = {1, 2};


Field[1] = Distance;
Field[1].CurvesList = {1};

Field[2] = MathEval;
Field[2].F = Sprintf("F1^1.25/12.+ %g", 0.2);
Background Field = 2;


Physical Curve("circle")  = {1};
Physical Curve("farfield") = {2};
Physical Surface("internal") = {1};


