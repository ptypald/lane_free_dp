/* 3 state variables -> x, y, vx */
/* controls -> long. acceleration and lateral speed */
#include <stdio.h>
#include <time.h>
#include <filesystem>
 
// include other headers
#include "dionysis_controller.h"
#include "data.h"
#include "controllers.h"
#include "vehicle.h"
#include "dddp.h"
#include "ddp.h"
#include "idm.h"

using namespace std;

void deleteDirectoryContents(const std::__fs::filesystem::path& dir) {
    for (const auto& entry : std::__fs::filesystem::directory_iterator(dir)) 
        std::__fs::filesystem::remove_all(entry.path()); 
}

/* Colors */
#define RED   "\x1B[1;31m"
#define GRN   "\x1B[1;32m"
#define RESET "\x1B[0m"

int main() {

	clock_t tic, toc;
	double cpu_time;
	
    // clear all output folders
    string path;
    // path = "../outputs/ddp/journal/solutions/ddp/temp"; deleteDirectoryContents(path);
    // path = "../outputs/feedback/solution"; deleteDirectoryContents(path);

	tic = clock();
	/* start procedure */

	// initialize and read data
	Data in;
	in.read();
    if (in.timeHorizon != 0.0)
        in.numsteps = in.timeHorizon / in.step;

	// Initialize vehicle(s)
	Vehicle veh(in);

    // use Dionysis controller for initial trajectory
    double d_step = 0.001;
    double tmax = in.numsteps*in.step;
    int d_numsteps = (int)tmax / d_step;
	Dionysis_Feedback fb(veh, d_numsteps, d_step, 1);
    fb.run(veh);

    // solve with DDP
    // -- Feedback initial
    DDP ddp(veh, in.numsteps, in.step, fb.opt_x, fb.opt_vx, fb.opt_y);
    int err = ddp.run(veh);
    if (err == 0)
        ddp.exportSimData(0);
    else if (err == 1)
        ddp.exportSimData(-1);

	/* end procedure */
	toc = clock();
	cpu_time = (double)(toc - tic) / CLOCKS_PER_SEC;
	fprintf(stderr, "-- Procedure Total CPU Time used: %.6f --\n", cpu_time);

	return 0;
}
