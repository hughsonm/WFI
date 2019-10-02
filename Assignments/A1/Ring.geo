ff = 1e9;
cc = 2.9979e8;
lambda = cc/ff;

cl = lambda/10;
ring_inner =   0.25*lambda;
ring_width = 0.05*lambda;

BKG_TAG = 1000;
RING_TAG = 1001;
BOUNDARY_TAG = 2000;
//+
Point(1) = {ring_inner, 0, 0, cl};
//+
Point(2) = {ring_inner+ring_width, 0, 0, cl};
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{1};
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{3};
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{5};
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{6};
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{2};
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{7};
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{8};
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{9};
}
//+
Curve Loop(1) = {6, 7, 8, 5};
//+
Curve Loop(2) = {2, 3, 4, 1};
//+
Plane Surface(1) = {1, 2};
//+
Physical Surface(1001) = {1};
