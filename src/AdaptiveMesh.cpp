#include "AdaptiveMesh.h"

AdaptiveMesh::AdaptiveMesh(int resolution) : resolution(resolution) {}

void AdaptiveMesh::refineAround(const Point& p) {
    // Placeholder implementation - just store the point
    meshNodes.push_back(p);
}
