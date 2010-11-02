////////////////////////////////////////////////////////////////////////////////
//
// Astronomical image alignment software. Copyright Edward Rosten 2009.
//

#include <gsl/gsl_randist.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <TooN/TooN.h>
#include <TooN/helpers.h>
#include <cvd/image_io.h>
#include <gvars3/instances.h>

using namespace std;
using namespace CVD;
using namespace TooN;
using namespace GVars3;


int main(int argc, char ** argv)
{
	GUI.LoadFile("ncc.cfg");
	GUI.parseArguments(argc, argv);

	ImageRef size = GV3::get<ImageRef>("sim.size", "[256 256]", -1);
	int nstars = GV3::get<int>("sim.num_stars", 30, -1);
	double max_offset = GV3::get<double>("sim.max_offset", 20, -1);
	double fwhm = GV3::get<double>("sim.star_fwhm", 30, -1);
	double snr = GV3::get<double>("sim.star_snr", 30, -1);
	double background = GV3::get<double>("sim.background", 30, -1);

	gsl_rng * eng = gsl_rng_alloc(gsl_rng_mt19937);

	double sigma2 = fwhm *fwhm/ (4 * log(2));

	//Scatter some stars about.
	vector<Vector<2> > stars;
	for(int i=0; i < nstars; i++)
		stars.push_back(makeVector(gsl_ran_flat(eng,0,size.x), gsl_ran_flat(eng,0,size.y)));


	Image<float> sim(size, 0);

	for(int y=0; y < sim.size().y; y++)
		for(int x=0; x < sim.size().x; x++)
		{
			//Find the contribution from each star
			double pix=1;
			Vector<2> p = makeVector(x,y);
			for(unsigned int i=0; i < stars.size(); i++)
				pix += exp(-norm_sq(p - stars[i])/(2*sigma2));

			sim[y][x] = gsl_ran_poisson(eng, background*pix);
		}

	img_save(sim, "foo.tiff");
}
