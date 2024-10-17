#include <iostream>
#include <cmath>
#include <vector>
#include <memory>
#include <fstream>
#include <algorithm>

// Define M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Constants
const double EPSILON = 1e-6;
const int MAX_DEPTH = 5;
const int IMAGE_WIDTH = 800;
const int IMAGE_HEIGHT = 600;

// Utility Functions
inline double degrees_to_radians(double degrees) {
    return degrees * M_PI / 180.0;
}

// Forward declaration for non-member operator*
struct Vector3;

// Non-member operator* for double * Vector3
inline Vector3 operator*(double scalar, const Vector3& v);

// Vector3 Class for 3D Vector Operations
struct Vector3 {
    double x, y, z;

    Vector3() : x(0), y(0), z(0){}
    Vector3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    // Vector addition
    Vector3 operator+(const Vector3& v) const { return Vector3(x + v.x, y + v.y, z + v.z); }

    // Vector subtraction
    Vector3 operator-(const Vector3& v) const { return Vector3(x - v.x, y - v.y, z - v.z); }

    // Unary minus
    Vector3 operator-() const { return Vector3(-x, -y, -z); }

    // Vector multiplication (element-wise)
    Vector3 operator*(const Vector3& v) const { return Vector3(x * v.x, y * v.y, z * v.z); }

    // Scalar multiplication
    Vector3 operator*(double scalar) const { return Vector3(x * scalar, y * scalar, z * scalar); }

    // Scalar division
    Vector3 operator/(double scalar) const { return Vector3(x / scalar, y / scalar, z / scalar); }

    // Dot product
    double dot(const Vector3& v) const { return x * v.x + y * v.y + z * v.z; }

    // Cross product
    Vector3 cross(const Vector3& v) const { 
        return Vector3(
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        );
    }

    // Magnitude of the vector
    double length() const { return std::sqrt(x*x + y*y + z*z); }

    // Normalize the vector
    Vector3 normalized() const { 
        double len = length();
        if (len == 0) return Vector3(0,0,0);
        return (*this) / len; 
    }
};

// Implementation of non-member operator* (double * Vector3)
inline Vector3 operator*(double scalar, const Vector3& v) { 
    return Vector3(v.x * scalar, v.y * scalar, v.z * scalar); 
}

// Ray Class
struct Ray {
    Vector3 origin;
    Vector3 direction;

    Ray(const Vector3& o, const Vector3& d) : origin(o), direction(d.normalized()) {}
};

// Color Class
struct Color {
    double r, g, b;

    Color() : r(0), g(0), b(0){}
    Color(double r_, double g_, double b_) : r(r_), g(g_), b(b_) {}

    // Clamp color values to [0,1]
    void clamp() {
        r = std::min(1.0, std::max(0.0, r));
        g = std::min(1.0, std::max(0.0, g));
        b = std::min(1.0, std::max(0.0, b));
    }

    // Multiply by scalar
    Color operator*(double scalar) const { return Color(r * scalar, g * scalar, b * scalar); }

    // Add two colors
    Color operator+(const Color& c) const { return Color(r + c.r, g + c.g, b + c.b); }

    // Element-wise multiplication
    Color operator*(const Color& c) const { return Color(r * c.r, g * c.g, b * c.b); }
};

// Material Class
struct Material {
    Color color;
    double reflectivity; // 0 = diffuse, 1 = perfect mirror

    Material(const Color& c = Color(1,1,1), double refl = 0.0) : color(c), reflectivity(refl) {}
};

// Abstract Object Class
struct Object {
    Material material;

    Object(const Material& m) : material(m) {}
    virtual bool intersect(const Ray& ray, double& t, Vector3& normal) const = 0;
};

// Sphere Class
struct Sphere : public Object {
    Vector3 center;
    double radius;

    Sphere(const Vector3& c, double r, const Material& m) : Object(m), center(c), radius(r) {}

    bool intersect(const Ray& ray, double& t, Vector3& normal) const override {
        Vector3 oc = ray.origin - center;
        double a = ray.direction.dot(ray.direction);
        double b = 2.0 * oc.dot(ray.direction);
        double c_val = oc.dot(oc) - radius * radius;
        double discriminant = b*b - 4*a*c_val;

        if (discriminant < 0) return false;
        double sqrt_disc = std::sqrt(discriminant);
        double t0 = (-b - sqrt_disc) / (2*a);
        double t1 = (-b + sqrt_disc) / (2*a);

        t = t0;
        if (t < EPSILON) {
            t = t1;
            if (t < EPSILON) return false;
        }

        Vector3 hit_point = ray.origin + ray.direction * t;
        normal = (hit_point - center).normalized();
        return true;
    }
};

// Plane Class
struct Plane : public Object {
    Vector3 point;   // A point on the plane
    Vector3 normal_vec; // Normal to the plane

    Plane(const Vector3& p, const Vector3& n, const Material& m) : Object(m), point(p), normal_vec(n.normalized()) {}

    bool intersect(const Ray& ray, double& t, Vector3& normal) const override {
        double denom = normal_vec.dot(ray.direction);
        if (std::abs(denom) < EPSILON) return false; // Parallel

        t = (point - ray.origin).dot(normal_vec) / denom;
        if (t < EPSILON) return false;

        normal = normal_vec;
        return true;
    }
};

// Light Class
struct Light {
    Vector3 position;
    Color color;
    double intensity;

    Light(const Vector3& pos, const Color& c, double inten) : position(pos), color(c), intensity(inten) {}
};

// Scene Class
struct Scene {
    std::vector<std::shared_ptr<Object>> objects;
    std::vector<Light> lights;
    Color background_color;

    Scene(const Color& bg = Color(0.2, 0.7, 1.0)) : background_color(bg) {}
};

// Function to compute lighting at a point
Color compute_lighting(const Scene& scene, const Vector3& point, const Vector3& normal, const Vector3& view_dir) {
    Color result(0, 0, 0);
    for (const auto& light : scene.lights) {
        Vector3 light_dir = (light.position - point).normalized();

        // Shadow check
        Ray shadow_ray(point + normal * EPSILON, light_dir);
        bool in_shadow = false;
        for (const auto& obj : scene.objects) {
            double t;
            Vector3 n;
            if (obj->intersect(shadow_ray, t, n)) {
                in_shadow = true;
                break;
            }
        }
        if (in_shadow) continue;

        // Diffuse shading
        double diff = std::max(0.0, normal.dot(light_dir));
        Color diffuse = light.color * diff * light.intensity;

        // Specular shading
        Vector3 reflect_dir = (2.0 * normal.dot(light_dir) * normal - light_dir).normalized();
        double spec = std::pow(std::max(0.0, reflect_dir.dot(view_dir)), 32);
        Color specular = light.color * spec * light.intensity;

        result = result + diffuse + specular;
    }
    return result;
}

// Trace Ray Function
Color trace_ray(const Scene& scene, const Ray& ray, int depth) {
    if (depth > MAX_DEPTH) return scene.background_color;

    double closest_t = INFINITY;
    std::shared_ptr<Object> hit_object = nullptr;
    Vector3 hit_normal;

    // Find closest intersection
    for (const auto& obj : scene.objects) {
        double t;
        Vector3 normal;
        if (obj->intersect(ray, t, normal)) {
            if (t < closest_t) {
                closest_t = t;
                hit_object = obj;
                hit_normal = normal;
            }
        }
    }

    if (!hit_object) return scene.background_color;

    Vector3 hit_point = ray.origin + ray.direction * closest_t;
    Vector3 view_dir = -ray.direction;

    // Compute local color
    Color local_color = hit_object->material.color * compute_lighting(scene, hit_point, hit_normal, view_dir);

    // Reflection
    if (hit_object->material.reflectivity > 0) {
        Vector3 reflect_dir = (ray.direction - hit_normal * 2.0 * ray.direction.dot(hit_normal)).normalized();
        Ray reflect_ray(hit_point + hit_normal * EPSILON, reflect_dir);
        Color reflected_color = trace_ray(scene, reflect_ray, depth + 1);
        local_color = local_color * (1 - hit_object->material.reflectivity) + reflected_color * hit_object->material.reflectivity;
    }

    local_color.clamp();
    return local_color;
}

// Function to write the image to a PPM file
void write_ppm(const std::string& filename, const std::vector<Color>& framebuffer, int width, int height) {
    std::ofstream ofs(filename, std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (const auto& color : framebuffer) {
        unsigned char r = static_cast<unsigned char>(std::min(1.0, color.r) * 255);
        unsigned char g = static_cast<unsigned char>(std::min(1.0, color.g) * 255);
        unsigned char b = static_cast<unsigned char>(std::min(1.0, color.b) * 255);
        ofs << r << g << b;
    }
    ofs.close();
}

// Main Function
int main() {
    // Initialize Scene
    Scene scene(Color(0.2, 0.7, 1.0)); // Sky blue background

    // Add Objects
    // Sphere: center, radius, material(color, reflectivity)
    scene.objects.push_back(std::make_shared<Sphere>(Vector3(0, 0, -5), 1.0, Material(Color(1, 0, 0), 0.5))); // Red sphere
    scene.objects.push_back(std::make_shared<Sphere>(Vector3(2, 0, -6), 1.0, Material(Color(0, 1, 0), 0.3))); // Green sphere
    scene.objects.push_back(std::make_shared<Sphere>(Vector3(-2, 0, -6), 1.0, Material(Color(0, 0, 1), 0.3))); // Blue sphere

    // Plane: point, normal, material
    scene.objects.push_back(std::make_shared<Plane>(Vector3(0, -1, 0), Vector3(0, 1, 0), Material(Color(1, 1, 1), 0.2))); // Ground

    // Add Lights
    scene.lights.emplace_back(Vector3(5, 5, 0), Color(1,1,1), 1.0);
    scene.lights.emplace_back(Vector3(-5, 5, 5), Color(1,1,1), 0.5);

    // Camera Setup
    Vector3 camera_pos(0, 0, 0);
    double fov = 90;
    double aspect_ratio = static_cast<double>(IMAGE_WIDTH) / IMAGE_HEIGHT;
    double scale = std::tan(degrees_to_radians(fov * 0.5));

    // Framebuffer
    std::vector<Color> framebuffer(IMAGE_WIDTH * IMAGE_HEIGHT);

    // Render Loop
    for (int j = 0; j < IMAGE_HEIGHT; ++j) {
        for (int i = 0; i < IMAGE_WIDTH; ++i) {
            // Normalize pixel coordinates to [-1,1]
            double x = (2 * (i + 0.5) / static_cast<double>(IMAGE_WIDTH) - 1) * aspect_ratio * scale;
            double y = (1 - 2 * (j + 0.5) / static_cast<double>(IMAGE_HEIGHT)) * scale;

            Vector3 dir(x, y, -1); // Assuming camera looks towards -Z
            Ray ray(camera_pos, dir);

            Color color = trace_ray(scene, ray, 0);
            framebuffer[j * IMAGE_WIDTH + i] = color;
        }
        // Optional: Print progress
        if (j % 50 == 0) std::cout << "Rendering row " << j << " of " << IMAGE_HEIGHT << std::endl;
    }

    // Write the framebuffer to a PPM file
    write_ppm("output.ppm", framebuffer, IMAGE_WIDTH, IMAGE_HEIGHT);
    std::cout << "Rendering completed. Image saved to 'output.ppm'." << std::endl;

    return 0;
}
