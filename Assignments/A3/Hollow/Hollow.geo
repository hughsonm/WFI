ff = 1e9;
cc = 3e8;
lambda = cc/ff;
cl = lambda / 10;
TARGET_TAG = 1000;

obj_rad = 2.1*lambda;
wall_thickness = 0.3*lambda;

dx = lambda/2.0;
dy = lambda/3.0;
//+
//+
Point(1) = {obj_rad, 0, 0, cl};
//+
Extrude {wall_thickness, 0, 0} {
  Point{1};
}
//+
Extrude {0, obj_rad, 0} {
  Curve{1};
}
//+
Extrude {{0, 0, 1}, {obj_rad, obj_rad, 0}, Pi/2} {
  Curve{2};
}
//+
Extrude {{0, 0, 1}, {obj_rad, obj_rad/2, 0}, Pi/2} {
  Curve{6};
}
//+
Extrude {{0, 0, 1}, {obj_rad/2, obj_rad/2, 0}, +Pi/2} {
  Curve{9};
}
//+
Extrude {{0, 0, 1}, {obj_rad/2, obj_rad/2-wall_thickness, 0}, -Pi/2} {
  Curve{13};
}
//+
Extrude {0, -(obj_rad/2-wall_thickness), 0} {
  Curve{16};
}
//+
Extrude {0, -obj_rad, 0} {
  Curve{1};
}
//+
Line(26) = {17, 18};
//+
Line(27) = {16, 18};
//+
Curve Loop(1) = {27, -26, -19};
//+
Plane Surface(27) = {1};
//+
Physical Surface(TARGET_TAG) = {27, 26, 5, 8, 12, 15, 18, 22};
