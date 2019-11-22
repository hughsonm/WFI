ff = 1e9;
cc = 2.9979e8;
wavelength = cc/ff;

cl = wavelength/10;

dom_rad = wavelength*2;
//+
Point(1) = {dom_rad, dom_rad, 0, cl};
//+
Extrude {-2*dom_rad, 0, 0} {
  Point{1};
}
//+
Extrude {0, -2*dom_rad, 0} {
  Curve{1};
}
//+
Physical Surface(1000) = {5};
