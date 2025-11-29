// Gmsh project created on Mon Nov 24 22:55:38 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {-0.05, 0, 0, 1.0};
//+
Point(2) = {0.05, 0, 0, 1.0};
//+
Point(3) = {-0.05, 0.01, 0, 1.0};
//+
Point(4) = {0.05, 0.01, 0, 1.0};
//+
Point(10) = {-0.04, 0.0, 0, 1.0};


//+
Line(1) = {1, 10};
//+
Line(2) = {10, 2};
//+
Line(3) = {2, 4};
//+
Line(4) = {3, 1};
//+
Point(5) = {-0.04, 0.0025, 0, 1.0};
//+
Point(6) = {-0.04, 0.0075, 0, 1.0};
//+
Point(7) = {-0.03, 0.0025, 0, 1.0};
//+
Point(8) = {-0.03, 0.0075, 0, 1.0};

//+
Line(5) = {6, 5};
//+
Line(6) = {5, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 6};
//+
Point(9) = {-0.04, 0.01, 0, 1.0};
//+
Line(9) = {9, 6};
//+
Physical Curve("wall", 9) = {1,10,11,4};
//+
Physical Curve("anode", 10) = {2,3};
//+
Physical Curve("cathode", 11) = {12,9, 8, 7, 6, 5};
// We first define a Distance field (`Field[1]') on
// Curves 2 and 4. This field returns the distance to
// (100 equidistant points on) curves 1 and 3.
Field[1] = Distance;
Field[1].CurvesList = {2, 18,19,20,21,3,4, 5, 6, 7, 9, 11,12, 17, 8, 10};
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
Field[2].SizeMin = 0.0003;
Field[2].SizeMax = 0.001;
Field[2].DistMin = 0.001;
Field[2].DistMax = 0.002;

// Then, we select Field 2 as the field determining the desired size of the mesh.
Background Field = 2;

// When the element size is fully specified by a mesh size field (as it is in
// this example), it is often desirable to set

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
//+
Line(10) = {4, 9};
//+
Line(11) = {9, 3};//+
Line(12) = {5, 10};
//+
Curve Loop(1) = {10, 9, -8, -7, -6, 12, 2, 3};
//+
Plane Surface(1) = {1};
