cl = 10e-2;
rad = 20e-2;

OMEGA_TAG = 1000;
BOUNDARY_TAG = 2000;
//+
Point(1) = {rad, 0, 0, cl};
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
Physical Surface(OMEGA_TAG) = {1};
//+
Physical Curve(BOUNDARY_TAG) = {2, 3, 4, 1};
