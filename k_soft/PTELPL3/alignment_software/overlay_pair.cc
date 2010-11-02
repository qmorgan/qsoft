////////////////////////////////////////////////////////////////////////////////
//
// Astronomical image alignment software. Copyright Edward Rosten 2009.
//

#include <cvd/image_io.h>
#include <cvd/cpu_hacks.h>
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
	enableFPE();
	GUI.LoadFile("mosaic.cfg");
	int lastarg = GUI.parseArguments(argc, argv);

	float centre=GV3::get<float>("centre", 0,-1);
	float range=GV3::get<float>("range", 0,-1);
	float zoom=GV3::get<float>("zoom", 3, 1);
	float border=GV3::get<float>("border", .2, 1);
	
	//Load all the unregistered images
	vector<string> image_names;
	for(int i=lastarg; i < argc; i++)
		image_names.push_back(argv[i]);

	if(image_names.size() != 2)
	{
		cerr << "Specify exactly 2 images\n";
		return 1;
	}

	Vector<2> mean;

	cin >> mean;

	Image<float> images = img_load(image_names[1]);
	Image<float> source = img_load(image_names[0]);

	ImageRef out_size = ir_rounded(vec(source.size()) * (1 + 2 * border) * zoom);

	Image<Rgb<byte> > out(out_size, Rgb<byte>(0,0,255));
	ImageRef out_off=ir_rounded(vec(source.size()) * border);

	//Convert and copy the original
	Image<byte> bytebit = convert_image(scale_image(source, centre, range));
	SubImage<Rgb<byte> >::iterator oi = out.sub_image(out_off, source.size()).begin();

	image_interpolate<Interpolate::Bicubic, byte> si(bytebit);

	for(int r=0; r < out.size().y; r++)
		for(int c=0; c < out.size().x; c++)
		{
			Vector<2> p = makeVector(c, r)/zoom - vec(out_off);
			if(si.in_image(p))
			{
				out[r][c].red = si[p];
				out[r][c].blue = 0;
			}
		}

	
	//Convert and copy the registered chip
	bytebit = convert_image(scale_image(images,centre,range));

	image_interpolate<Interpolate::Bicubic, byte> bi(bytebit);

	for(int r=0; r < out.size().y; r++)
		for(int c=0; c < out.size().x; c++)
		{
			Vector<2> p = makeVector(c, r)/zoom - vec(out_off) - mean;
			if(bi.in_image(p))
			{
				out[r][c].green = bi[p];
				out[r][c].blue = 0;
			}
		}

	//Draw on the box match positions
	img_save(out, cout, ImageType::PNG);

}
catch(Exceptions::All a)
{
	cerr << "Error: " << a.what << endl;
}
}
