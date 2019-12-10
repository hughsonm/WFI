ff = 1e9;
cc = 3e8;
wavelength = cc/ff;

cl = wavelength/10;

target_rad = 0.5*wavelength;

x_offset = 0.1*wavelength;
y_offset = -0.2*wavelength;
//+
Point(1) = {x_offset+target_rad, y_offset, 0, cl};
//+
Rotate {{0, 0, 1}, {x_offset, y_offset, 0}, 2*Pi/3} {
  Point{1};
}
//+
Point(2) = {x_offset+target_rad, y_offset, 0, cl};
//+
Rotate {{0, 0, 1}, {x_offset, y_offset, 0}, -2*Pi/3} {
  Point{2};
}
//+
Point(3) = {x_offset+target_rad, y_offset, 0, cl};
//+
Line(1) = {3, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Curve Loop(1) = {2, 3, 1};
//+
Plane Surface(1) = {1};
//+
Physical Surface(1000) = {1};
