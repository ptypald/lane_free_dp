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
    path = "../outputs/dddp/domain"; deleteDirectoryContents(path);
    path = "../outputs/dddp/solution"; deleteDirectoryContents(path);
    path = "../outputs/ddp/solution"; deleteDirectoryContents(path);
    path = "../outputs/idm/solution"; deleteDirectoryContents(path);
    path = "../outputs/feedback/solution"; deleteDirectoryContents(path);

	tic = clock();
	/* start procedure */

	// initialize and read data
	Data in;
	in.read();
    if (in.timeHorizon != 0.0)
        in.numsteps = in.timeHorizon / in.step;

	// Initialize vehicle(s)
	Vehicle veh(in);

    // IDM
    // IDM idm(veh, in.numsteps, in.step);
    // idm.idm_run();
    // idm.exportSimData(0);

    // use Dionysis controller for initial trajectory
    double d_step = 0.001;
    double tmax = in.numsteps*in.step;
    int d_numsteps = (int)tmax / d_step;
	Dionysis_Feedback fb(veh, d_numsteps, d_step, 1);
    fb.run(veh);
    fb.exportSimData(0);

    // solve with DDDP
    DDDP dddp(veh, in.numsteps, in.step, fb.opt_x, fb.opt_vx, fb.opt_y);
    // dddp.run(veh);
    // TODO: change discretizexv() --> use fmod? to get the closer grid point, in case interpolation is needed (zero interpolation) [for T != 1.0]
    dddp.exportSimData(-1);
    // return 0;
    
    // solve with DDP
    // -- IDM initial
    // DDP ddp(veh, in.numsteps, in.step, idm.opt_x, idm.opt_vx, idm.opt_y, idm.opt_ux, idm.opt_uy);    
    // -- Feedback initial
    DDP ddp(veh, in.numsteps, in.step, dddp.init_x, dddp.init_vx, dddp.init_y, dddp.init_ux, dddp.init_uy);    
    // -- DDP initial
    // DDP ddp(veh, in.numsteps, in.step, dddp.opt_x, dddp.opt_vx, dddp.opt_y, dddp.opt_ux, dddp.opt_uy);
    // DDP ddp(veh, in.numsteps, in.step, testInit_x, testInit_vx, testInit_y, testInit_vy, testInit_ux, testInit_uy);
    ddp.run(veh);
    ddp.exportSimData(0);

	/* end procedure */
	toc = clock();
	cpu_time = (double)(toc - tic) / CLOCKS_PER_SEC;
	fprintf(stderr, "-- Procedure Total CPU Time used: %.6f --\n", cpu_time);

	return 0;
}
