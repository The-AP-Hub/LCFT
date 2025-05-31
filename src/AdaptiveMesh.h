#ifndef LCFT_ADAPTIVE_MESH_H
#define LCFT_ADAPTIVE_MESH_H

#include "Point.h"
#include <vector>

class AdaptiveMesh {
public:
    AdaptiveMesh(int resolution);
    void refineAround(const Point& p);
    const std::vector<Point>& nodes() const { return meshNodes; }
private:
    int resolution;
    std::vector<Point> meshNodes; // placeholder for mesh nodes
};

#endif // LCFT_ADAPTIVE_MESH_H
