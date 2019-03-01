#include <chrono>

#include "Model.h"
//#include "Burgers.h"

using namespace std;

int main(int argc, char* argv[]) {
    Model m(argc, argv);
    if (m.isVerbose()) m.printParameters();
    if (m.isValid()) {
        cout << "Now running the model with given parameters." << endl;
//                Burgers b(m);

        // Call code to initialise the problem here

        typedef std::chrono::high_resolution_clock hrc;
        typedef std::chrono::milliseconds ms;
        hrc::time_point start = hrc::now();

        // Call code to perform time integration here

        hrc::time_point end = hrc::now();

        // Calculate final energy and write output

    } else {
        cerr << "The model was not run as it is not valid with the given parameters." << endl
             << "Please choose a different set of parameters and try again or use '-h' flag for the help message."
             << endl;
    }

    return 0;
}