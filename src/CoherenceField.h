#ifndef LCFT_COHERENCE_FIELD_H
#define LCFT_COHERENCE_FIELD_H

#include "Point.h"
#include <vector>

class CoherenceField {
public:
    CoherenceField();
    double valueAt(const Point& p) const;
private:
    std::vector<double> fieldValues; // placeholder for field data
};

#endif // LCFT_COHERENCE_FIELD_H
