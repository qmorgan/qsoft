#include <cvd/image_io.h>
#include <vector>
#include <iterator>
#include <gvars3/instances.h>
using namespace CVD;
using namespace std;
using namespace GVars3;


Image<vector<ImageRef> > missing_pixels_interal_images(const Image<bool>& goodpixels)
{

}

int main(int argc, char** argv)
{
	GUI.parseArguments(argc, argv);

	Image<bool> goodpixels = img_load(cin);

	Image<vector<ImageRef> > region(goodpixels.size());

	vector<ImageRef> r;
	//Do the bottom row
	for(int x=0; x < region.size().x; x++)
	{
		if(!goodpixels[0][x])
			r.push_back(ImageRef(x,0));
		region[0][x] = goodpixels;
	}

	for(int y=1; y < region.size().y; y++)
	{
		r.clear();

		for(int x=0; x < region.size().x; x++)
		{
			if(!goodpixels[y][x])
				r.push_back(ImageRef(x,y));

			region[y][x] = region[y-1][x];
			copy(r.begin(), r.end(), back_inserter(region[y][x]));
		}
	}

}

