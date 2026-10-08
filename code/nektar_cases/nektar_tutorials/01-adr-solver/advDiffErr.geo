// Bigger square points
Point(1) = {-1.00, -1.00,   0.00, 0.4};
Point(2) = { 1.00, -1.00,   0.00, 0.4};
Point(3) = { 1.00,  1.00,   0.00, 0.4};
Point(4) = {-1.00,  1.00,   0.00, 0.4};

// Lines of the square
Line(1)  = {1, 2};
Line(2)  = {2, 3};
Line(3)  = {3, 4};
Line(4)  = {4, 1};

//Line Loops to create physical surfaces
Line Loop(1)        = {1, 2, 3, 4};
Plane Surface(1)    = {1};
Physical Surface(1) = {1};


Transfinite Line {1, -3} = 11 Using Progression 1;
Transfinite Line {2, -4} = 11 Using Progression 1;
Transfinite Surface{1};
Recombine Surface{1};


// Defining physical lines for applying the boundary
// conditions
Physical Line(100) = {1};     // botoom
Physical Line(200) = {2};     // right
Physical Line(300) = {3};     // top
Physical Line(400) = {4};     // left
