#ifndef BURGERS_H
#define BURGERS_H

#include <iostream>
#include <string>
#include "Model.h"
#include <blaze/Math.h>
using namespace std;


class Burgers {
    public:
        explicit Burgers(Model *m_);
        void initializeVelocityField();


    private:
        Model* m;
        blaze::DynamicMatrix<double,blaze::rowMajor> u, v;
        double x(size_t col), y(size_t row);
        double r_thresh = 1;
};

#endif