# Ray Tracer

## Overview

**Ray Tracer** is a single-file C++ application that implements a basic ray tracing algorithm to render 3D scenes. This project demonstrates advanced C++ programming techniques, including object-oriented design, recursion, and optimization strategies. The ray tracer supports multiple geometric objects, lighting models, shadows, and reflections, producing realistic rendered images.

## Features

- **Geometric Objects:**
    
    - Spheres
    - Planes
- **Lighting:**
    
    - Point light sources with configurable position, color, and intensity
    - Diffuse and specular shading
    - Shadows
- **Materials:**
    
    - Diffuse and reflective surfaces with configurable reflectivity
- **Camera:**
    
    - Configurable position and field of view
- **Rendering:**
    
    - Recursive ray tracing with reflection depth control
    - Image output in PPM format (can be converted to PNG or JPEG)

### Viewing the Output

The PPM format is widely supported, but you can convert it to more common formats like PNG using tools like ImageMagick:

## Customization

- **Scene Configuration:**
    
    - Modify the `main` function in `ray_tracer.cpp` to add or remove objects and lights.
    - Adjust material properties like color and reflectivity to change object appearances.
- **Rendering Settings:**
    
    - Change `IMAGE_WIDTH` and `IMAGE_HEIGHT` constants to adjust the resolution.
    - Modify `MAX_DEPTH` to control reflection recursion depth.


## Acknowledgments

- Inspired by fundamental ray tracing principles in computer graphics.
- Utilizes C++11 features for efficient and modern code design.
