#include "SystemMatrices.h"

SystemMatrices::SystemMatrices(int size) : size(size), data(size*size, 0.0) {}

double SystemMatrices::getValue(int i, int j) const {
    return data[i*size + j];
}

void SystemMatrices::setValue(int i, int j, double value) {
    data[i*size + j] = value;
}
