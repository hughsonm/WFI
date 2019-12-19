ff = 1e9;
cc = 3e8;
wavelength = cc/ff;

cl = wavelength/10;

target_rad = 0.5*wavelength;

x_offset = 0.1*wavelength;
y_offset = -0.2*wavelength;

wall_width = 0.1*wavelength;
ball_rad = 0.15*wavelength;
//+

//+
//+
Point(2) = {target_rad+wall_width, target_rad+wall_width, 0, cl};
//+
Point(3) = {-(target_rad+wall_width), target_rad+wall_width, 0, cl};
//+
Point(4) = {-(target_rad+wall_width), -(target_rad+wall_width), 0, cl};
//+
Point(5) = {(target_rad+wall_width), -(target_rad+wall_width), 0, cl};
//+
Point(6) = {(target_rad), (target_rad), 0, cl};
//+
Point(7) = {-(target_rad), (target_rad), 0, cl};
//+
Point(8) = {-(target_rad), -(target_rad), 0, cl};
//+
Point(9) = {(target_rad), -(target_rad), 0, cl};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 2};
//+
Line(5) = {6, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 9};
//+
Line(8) = {9, 6};
//+
Point(10) = {x_offset+ball_rad, y_offset, 0, cl};
//+
Extrude {{0, 0, 1}, {x_offset, y_offset, 0}, Pi/3*2} {
  Point{10};
}
//+
Extrude {{0, 0, 1}, {x_offset, y_offset, 0}, Pi/3*2} {
  Point{11};
}
//+
Extrude {{0, 0, 1}, {x_offset, y_offset, 0}, Pi/3*2} {
  Point{13};
}
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(1) = {1, 2};
//+
Curve Loop(3) = {9, 10, 11};
//+
Plane Surface(2) = {3};
//+
Physical Surface(1000) = {1};
//+
Physical Surface(1001) = {2};
