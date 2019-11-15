ff = 1e9;
cc = 3e8;
lambda = cc/ff;
cl = lambda / 10;
TARGET_TAG = 1000;

obj_rad = 1.2*lambda;
wall_thickness = 0.2*lambda;

dx = lambda/2.0;
dy = lambda/3.0;
//+
Point(1) = {obj_rad, 0, 0, cl};
//+
Extrude {wall_thickness, 0, 0} {
  Point{1};
}
//+
Extrude {0, 3*obj_rad, 0} {
  Curve{1};
}
//+
Extrude {{0, 0, 1}, {0, 3*obj_rad, 0}, Pi/2} {
  Curve{2};
}
//+
Extrude {{0, 0, 1}, {0, 3*obj_rad, 0}, Pi/2} {
  Curve{6};
}
//+
Extrude {0, -3*obj_rad, 0} {
  Curve{10};
}
//+
Extrude {{0, 0, 1}, {-obj_rad-wall_thickness, 0, 0}, -Pi/2} {
  Curve{14};
}
//+
Extrude {{0, 0, 1}, {-obj_rad-wall_thickness, -obj_rad-wall_thickness, 0}, Pi/2} {
  Curve{18};
}
//+
Extrude {{0, 0, 1}, {-obj_rad-wall_thickness, -obj_rad-wall_thickness, 0}, Pi/2} {
  Curve{21};
}
//+
Extrude {{0, 0, 1}, {-obj_rad-wall_thickness, -obj_rad-wall_thickness, 0}, Pi/4} {
  Curve{25};
}
//+
//+
Extrude {{0, 0, 1}, {obj_rad+wall_thickness, 0, 0}, Pi/2} {
  Curve{1};
}
//+
Extrude {{0, 0, 1}, {obj_rad+wall_thickness, -obj_rad-wall_thickness, 0}, -Pi/2} {
  Curve{33};
}
//+
Extrude {{0, 0, 1}, {obj_rad+wall_thickness, -obj_rad-wall_thickness, 0}, -Pi/2} {
  Curve{36};
}
//+
Extrude {{0, 0, 1}, {obj_rad+wall_thickness, -obj_rad-wall_thickness, 0}, -Pi/4} {
  Curve{40};
}
//+
Extrude {{0, 0, 1}, {0, -2*(obj_rad+wall_thickness), 0}, Pi/2} {
  Curve{44};
}
//+
Physical Surface(TARGET_TAG) = {51, 47, 43, 39, 35, 5, 9, 13, 17, 20, 24, 28, 32};
