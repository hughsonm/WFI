ff = 1e9;
cc = 3e8;
lambda = cc/ff;
cl = lambda / 10;

target_rad = 0.3*lambda;

TARGET_TAG = 1000;
//+
Point(1) = {target_rad, 0, 0, cl};
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{1};
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{2};
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{4};
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {
  Point{5};
}
//+
Curve Loop(1) = {2, 3, 4, 1};
//+
Plane Surface(1) = {1};
//+
Physical Surface(TARGET_TAG) = {1};
