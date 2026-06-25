
#include <iostream>
#include <string>

#include "MeshingVoroCrust.h"

int main(const int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    MeshingVoroCrustOptions options;
    options.input_mesh_files.push_back("OneBox.obj");
    options.r_max = 10000.0;
    options.Lip_const = 0.25;
    options.vc_ang_tol = 60.0;
    options.num_threads = 1;

    MeshingVoroCrust vc(options);
    int rc = vc.execute();

    if (rc != 0)
    {
        std::cerr << "ERROR: MeshingVoroCrust::execute() returned " << rc << std::endl;
        return 1;
    }

    std::cout << "MeshingVoroCrust options-based test passed." << std::endl;
    return 0;
}
