cc = 3e8;
ff = 1e9;
lambda = cc/ff;

cl = lambda/10;

cell_size = lambda/3;
//+
Point(1) = {0, 0, 0, cl};
//+
Extrude {-cell_size, 0, 0} {
  Point{1};
}
//+
Extrude {0, cell_size, 0} {
  Point{2};
}
//+
Extrude {2*cell_size, 0, 0} {
  Point{3};
}
//+
Extrude {0, cell_size, 0} {
  Point{4};
}
//+
Extrude {-3*cell_size, 0, 0} {
  Point{5};
}
//+
Extrude {0, -4*cell_size, 0} {
  Point{6};
}
//+
Extrude {cell_size, 0, 0} {
  Point{7};
}
//+
Extrude {0, cell_size, 0} {
  Point{8};
}
//+
Extrude {cell_size, 0, 0} {
  Point{9};
}
//+
Extrude {0, cell_size, 0} {
  Point{10};
}
//+
Curve Loop(1) = {5, 6, 7, 8, 9, 10, 1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Surface(1000) = {1};
