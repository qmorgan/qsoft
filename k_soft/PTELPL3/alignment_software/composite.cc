////////////////////////////////////////////////////////////////////////////////
//
// Astronomical image alignment software. Copyright Edward Rosten 2009.
//

#include <cvd/image_io.h>
#include <cvd/integral_image.h>
#include <cvd/draw.h>
#include <cvd/image_convert.h>
#include <cvd/vector_image_ref.h>
#include <cvd/image_interpolate.h>


#include <TooN/TooN.h>
#include <TooN/SymEigen.h>

#include <gvars3/instances.h>

#include <iostream>
#include <iterator>
#include <algorithm>
#include <set>
#include <cmath>

using namespace CVD;
using namespace GVars3;
using namespace std;
using namespace TooN;

float lim(float f)
{
	if(f < 0)
		return 0;
	else if(f > 1)
		return 1;
	else return f;
}

Image<float> scale_image(const Image<float>& i, float c, float r)
{
	Image<float> s(i.size());
	
	Image<float>::const_iterator a=i.begin();
	Image<float>::iterator b=s.begin();

	for(; a!= i.end(); a++, b++)
		*b = lim((*a-c)/r + 0.5);

	return s;
}

int main(int argc, char** argv)
{
try
{	
	GUI.LoadFile("mosaic.cfg");
	GUI.parseArguments(argc, argv);

	float zoom=GV3::get<float>("zoom", 3, 1);
	float border=GV3::get<float>("border", .2, 1);
	
	//Load all the unregistered images
	vector<string> image_names;
	vector<Vector<2> > offsets;

	for(;;)
	{
		Vector<2> o;
		string n;
		cin >> o >> n;

		if(cin.fail())
			break;
		
		image_names.push_back(n);
		offsets.push_back(o);

	}

	Image<vector<float> > median;
	ImageRef out_size;
	Vector<2> out_off=Zeros;

	for(unsigned int i=0; i < image_names.size(); i++)
	{
		Image<float> source = img_load(image_names[i]);
		Vector<2> mean = offsets[i];

		if(i == 0)
		{
			out_size = ir_rounded(vec(source.size()) * (1 + 2 * border) * zoom);
			out_off=vec(source.size()) * border;

			median.resize(out_size);
		}

		image_interpolate<Interpolate::Bicubic, float> si(source);

		for(int r=0; r < median.size().y; r++)
			for(int c=0; c < median.size().x; c++)
			{
				Vector<2> p = makeVector(c, r)/zoom - out_off - mean;
				if(si.in_image(p))
					median[r][c].push_back(si[p]);
			}

		cerr << "Composited " << i << " of " << image_names.size() << endl;
	}
	
	Image<float> out(median.size(), 0);

	for(int r=0; r < median.size().y; r++)
		for(int c=0; c < median.size().x; c++)
			if(!median[r][c].empty())
			{
				sort(median[r][c].begin(), median[r][c].end());
				out[r][c] = median[r][c][median[r][c].size()/2];
			}

	//Draw on the box match positions
	img_save(out, cout, ImageType::FITS);

}
catch(Exceptions::All a)
{
	cerr << "Error: " << a.what << endl;
}
}
