#ifndef BURGERS_H
#define BURGERS_H

#include <iostream>
#include <string>
#include "Model.h"
using namespace std;


class Burgers {
    public:
        Burgers(Model m);
        void initializeVelocityField();
};

#endif