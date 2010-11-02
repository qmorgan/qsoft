////////////////////////////////////////////////////////////////////////////////
//
// Astronomical image alignment software. Copyright Edward Rosten 2009.
//

#include <TooN/TooN.h>
#include <vector>
#include <gvars3/instances.h>

using namespace std;
using namespace TooN;
using namespace GVars3;

int main(int argc, char** argv)
{
	GUI.LoadFile("mosaic.cfg");
	GUI.parseArguments(argc, argv);

	double noise_sig_1 = GV3::get<double>("weighted_mean.noise_std_start", 0., -1);
	double noise_sig_2 = GV3::get<double>("weighted_mean.noise_std_end", 0., -1);
	int noise_var_steps = GV3::get<int>("weighted_mean.noise_steps", -1, -1);
	double tol = GV3::get<double>("weighted_mean.tolerance", 0., -1);

	Vector<2> mean = GV3::get<Vector<2> >("start", Zeros, -1);
	
	vector<Vector<2> > v;
	//Read in the data
	for(;;)
	{
		Vector<2> i;
		cin >> i;
		if(!cin.good())
			break;
		v.push_back(i);
	}

	if(v.size() == 0)
	{
		cout << "-1e99 -1e99\n";
		return 1;
	}

	double noise_std = noise_sig_1;
	double alpha = pow(noise_sig_1/noise_sig_2, 1./(noise_var_steps-1));
	
	for(int i=0; i < noise_var_steps; i++)
	{
		for(;;)
		{
			Vector<2> newmean = Zeros;
			double weights=0;

			for(unsigned int i=0; i < v.size(); i++)
			{
				double w = 1./(noise_std * noise_std + norm_sq(mean - v[i]));
				newmean +=  w * v[i];
				weights +=  w;
			}

			newmean /= weights;


			if(norm(newmean - mean) < tol)
				break;
			mean = newmean;
		}

		noise_std /= alpha;
	}

	cout << mean << endl;
}
