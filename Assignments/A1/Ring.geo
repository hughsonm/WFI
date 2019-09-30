ff = 1e9;
cc = 2.9979e8;
lambda = cc/ff;

cl = 10e-2;
bd_rad = 20e-2;
ring_rad =   0.275*lambda;
ring_width = 0.05*lambda;

BKG_TAG = 1000;
RING_TAG = 1001;
BOUNDARY_TAG = 2000;
//+
Point(1) = {bd_rad, 0, 0, cl};
//+
//+
Point(2) = {ring_rad-ring_width/2, 0, 0, cl};
//+
Point(3) = {ring_rad+ring_width/2, 0, 0, cl};
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{2}; Point{3}; Point{1};
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{7}; Point{6}; Point{4};
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{10}; Point{9}; Point{8};
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{11}; Point{12}; Point{13};
}
//+
Curve Loop(1) = {1, 6, 7, 10};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 5, 8, 11};
//+
Plane Surface(2) = {1, 2};
//+
Curve Loop(3) = {3, 4, 9, 12};
//+
Plane Surface(3) = {2, 3};
//+
Physical Surface(BKG_TAG) = {1, 3};
//+
Physical Surface(1001) = {2};
//+
Physical Curve(2000) = {3, 4, 9, 12};
