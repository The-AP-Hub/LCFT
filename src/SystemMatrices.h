#ifndef LCFT_SYSTEM_MATRICES_H
#define LCFT_SYSTEM_MATRICES_H

#include <vector>

class SystemMatrices {
public:
    SystemMatrices(int size);
    double getValue(int i, int j) const;
    void setValue(int i, int j, double value);
private:
    int size;
    std::vector<double> data; // flattened matrix
};

#endif // LCFT_SYSTEM_MATRICES_H
