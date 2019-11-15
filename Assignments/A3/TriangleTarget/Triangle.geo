ff = 1e9;
cc = 3e8;
lambda = cc/ff;
cl = lambda / 10;
TARGET_TAG = 1000;

tri_rad = 1.2*lambda;
dx = lambda/2.0;
dy = lambda/3.0;
//+
Point(1) = {tri_rad+dx, dy, 0, cl};
//+
Rotate {{0, 0, 1}, {dx, dy, 0}, Pi*2/3} {
  Point{1};
}
//+
Point(2) = {tri_rad+dx, dy, 0, cl};
//+
Rotate {{0, 0, 1}, {dx, dy, 0}, -Pi*2/3} {
  Point{2};
}
//+
Point(3) = {tri_rad+dx, dy, 0, cl};
//+
Line(1) = {3, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Curve Loop(1) = {1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Surface(TARGET_TAG) = {1};
