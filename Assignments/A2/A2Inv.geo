ff = 5e9;
cc = 3e8;

lambda = cc/ff;

cl = lambda / 10;

img_dom_rad = 7.2e-2;
IMG_TAG = 2000;
//+
Point(1) = {img_dom_rad, 0, 0, cl};
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
Physical Surface(IMG_TAG) = {1};
