//
// Created by raedler on 26.03.2026.
//

#ifndef GATE_MULTI_DETECTOR_POST_PROCESSING_VEC3_H
#define GATE_MULTI_DETECTOR_POST_PROCESSING_VEC3_H

// #include <cmath>

struct Vec3 {
    double x, y, z;
};

// subtraction
inline Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

// dot product
inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// cross product
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

// squared norm
inline double norm2(const Vec3& v) {
    return dot(v, v);
}

// norm
inline double norm(const Vec3& v) {
    return std::sqrt(norm2(v));
}

#endif //GATE_MULTI_DETECTOR_POST_PROCESSING_VEC3_H