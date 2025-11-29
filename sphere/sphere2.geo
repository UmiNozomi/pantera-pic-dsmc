// Gmsh project created on Fri Oct 24 23:57:17 2025
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.09, 0, Pi};
//+
Circle(2) = {0, 0.03, 0, 0.0025, 0, 2*Pi};
//+
Circle(3) = {0, 0.06, 0, 0.0025, 0, 2*Pi};
//+
Rotate {{0, 0, 1}, {0, 0, 0}, -Pi/5} {
  Duplicata { Curve{3}; Curve{2}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, -Pi/5} {
  Duplicata { Curve{4}; Curve{5}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/5} {
  Duplicata { Curve{3}; Curve{2}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/5} {
  Duplicata { Curve{8}; Curve{9}; }
}
//+
Line(12) = {2, 1};
//+
Curve Loop(1) = {12, 1};
// We first define a Distance field (`Field[1]') on
// Curves 2 and 4. This field returns the distance to
// (100 equidistant points on) curves 1 and 3.
Field[1] = Distance;
Field[1].CurvesList = {2, 18,19,20,21,3,4, 5, 6, 7, 9, 11, 17, 16,1, 8, 10,14,13,15};
Field[1].Sampling = 200;
// We then define a `Threshold' field, which uses the return value of the
// `Distance' field 1 in order to define a simple change in element size
// depending on the computed distances
//
// SizeMax -                     /------------------
//                              /
//                             /
//                            /
// SizeMin -o----------------/
//          |                |    |
//        Point         DistMin  DistMax
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.0008;
Field[2].SizeMax = 0.0015;
Field[2].DistMin = 0.0015;
Field[2].DistMax = 0.003;

// Then, we select Field 2 as the field determining the desired size of the mesh.
Background Field = 2;

// When the element size is fully specified by a mesh size field (as it is in
// this example), it is often desirable to set

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

//+
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/10} {
  Duplicata { Curve{3}; Curve{8}; }
}
//+
Rotate {{0, 0, 1}, {0, 0, 0}, -Pi/10} {
  Duplicata { Curve{3}; Curve{4}; }
}
//+
Curve Loop(2) = {1, 12};
//+
Curve Loop(3) = {10};
//+
Curve Loop(4) = {14};
//+
Curve Loop(5) = {8};
//+
Curve Loop(6) = {13};
//+
Curve Loop(7) = {3};
//+
Curve Loop(8) = {15};
//+
Curve Loop(9) = {4};
//+
Curve Loop(10) = {16};
//+
Curve Loop(11) = {6};
//+
Curve Loop(12) = {11};
//+
Curve Loop(13) = {9};
//+
Curve Loop(14) = {2};
//+
Curve Loop(15) = {5};
//+
Curve Loop(16) = {7};
//+
Plane Surface(1) = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
//+
Physical Curve("wall", 17) = {1};
//+
Physical Curve("anode", 18) = {10, 14, 8, 13, 3, 15, 4, 16, 6};
//+
Physical Curve("cathode", 19) = {11, 9, 2, 5, 7};
