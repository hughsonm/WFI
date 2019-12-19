ff = 1e9;
cc = 2.9979e8;
wavelength = cc/ff;

cl = wavelength/10;

dom_rad = 20e-2;
//+
Point(1) = {dom_rad, dom_rad, 0, cl*4};
Point(2) = {dom_rad, -dom_rad, 0, cl*2};
Point(3) = {-dom_rad, -dom_rad, 0, cl/2};
Point(4) = {-dom_rad, dom_rad, 0, cl/4};
//+
Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Surface(1000) = {1};
