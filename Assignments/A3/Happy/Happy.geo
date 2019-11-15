ff = 1e9;
cc = 3e8;
lambda = cc/ff;
cl = lambda / 10;
TARGET_TAG = 1000;

obj_rad = 0.5*lambda;
wall_thickness = 0.15*lambda;
//+
Point(1) = {obj_rad, 0, 0, cl};
//+
Point(2) = {-obj_rad, 0, 0, cl};
//+
Extrude {wall_thickness, 0, 0} {
  Point{1};
}
//+
Extrude {-wall_thickness, 0, 0} {
  Point{2};
}
//+
Extrude {0, 2*obj_rad, 0} {
  Curve{2}; Curve{1};
}
//+
Point(9) = {-obj_rad, -obj_rad, 0, cl};
//+
Extrude {-wall_thickness, 0, 0} {
  Point{9};
}
//+
Extrude {{0, 0, 1}, {0, -obj_rad, 0}, Pi/2} {
  Curve{11};
}
//+
Extrude {{0, 0, 1}, {0, -obj_rad, 0}, Pi/2} {
  Curve{12};
}
//+
Physical Surface(TARGET_TAG) = {6, 10, 19, 15};
