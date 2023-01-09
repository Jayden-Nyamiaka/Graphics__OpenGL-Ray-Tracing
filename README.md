# CS 171 Raytracer

## Controls

There are a number of controls you can use to do different things. Here's an overview of the availale controls:

- `t`: Toggle wireframe mode (useful for IOTest checking).
  - Renders shapes as wireframes rather than solid objects. Lighting won't look right with this on, but it's useful for IOTest.
- `o`: Toggle the IOTest.
  - Shows a number of small spheres which are rendered only if the IOTest function returns `true` at the given coordinate.
- `i`: Toggle the intersect test.
  - Renders a single ray and the normal at the intersection point as a pair of white lines.
- `c`: Toggle scene camera mode.
  - When enabled, the camera gets locked to the camera properties specified by the scene. There may be subtle differences between this camera and the actual render, but they should be mostly the same. Useful for checking your renders.
- `r`: Runs the raytrace.
  - Raytraces the scene. Utilizes a second thread, so your display should remain responsive. Do not to spam this.
- `q`: Exits the program.

## Code Overview

All the code you have to write for this assignment is contained in `assignment.cpp`, although it may be useful to look at other parts of the code from time to time, so here's a general overview of what the rest of the code does.

- `camera.h`/`camera.cpp`: Implements the Camera class.
- `glinclude.h`: Contains the imports needed for OpenGL to function.
- `image.h`/`image.cpp`: Implements PNG exporting (using libpng).
- `light.h`/`light.cpp`: Implements the Light class.
- `object.h`/`object.cpp`: Implements some utility code for the Object/Superquadric/Assembly classes.
- `parsing.h`: Code that parses the scene file (using libyaml-cpp).
- `renderer.cpp`: Main application code.
- `scene.h`/`scene.cpp`: Implements the Scene class.
- `transform.h`/`transform.cpp`: Implements the different transformations.
- `util.h`/`util.cpp`: Implements various utility functions (not useful to you).

For the most part, the only things you should care about are the `Camera`, `Light`, `Object`/`Superquadric`/`Assembly`, and `Scene` classes. There are overviews of what these different classes do below.

### `Camera`

A simple class that contains the information about the scene's camera.

### `Light`

A simple class that contains the information about a light in the scene.

### `Object`

A class that can represent either an `Assembly` or a `Superquadric`. Essentially, this class defines the interface that we want to have for both `Assembly` and `Superquadric`. All you should care about is the `IOTest` function, the `ClosestIntersection` function, and the list of transformations.

### `Superquadric`

A subclass of `Object` that represents a superquadric with a list of transformations applied to it. All you should care about is the `exp0`, `exp1`, and `mat` members, the functions mentioned above, and the `GetNormal` function, which returns the normal vector at a point on the surface of the superquadric.

### `Assembly`

A subclass of `Object` that represents a combination of several sub-objects (that can be either assemblies or superquadrics) with a list of transformations applied to it. All you care about here is the list of children and the functions mentioned in the `Object` class.

### `Scene`

A class that represents a loaded scene. You should pay attention to the member variables `camera` and `lights`, and the `ClosestIntersection` function, which returns the first intersection between the given ray and any of the objects in the scene (and relies on your implementation of `ClosestIntersection` on the `Assembly` and `Superquadric` classes).



- Added methods to object class that utlilizes transforms to get the transform matrices
put function signatures 
getForwardTransformMatrix
getInverseTransformMatrix

- Added methods to Superquadric:
// Helpful functions to compute the inside-outside function
double IOFunction(Vector3d pos);  
Vector3d IOGradient(Vector3d pos);

- need to do 1/20 of a PIXEL for close enough bound
- mention setting this manually on first run of RayTrace bc that comes before anything that needs it and pixel calcultions are done there
-remember to paste in pictures
