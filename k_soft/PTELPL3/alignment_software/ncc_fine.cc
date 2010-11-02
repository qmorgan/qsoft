////////////////////////////////////////////////////////////////////////////////
//
// Astronomical image alignment software. Copyright Edward Rosten 2009.
//

#include <cvd/image_io.h>
#include <cvd/interpolate.h>
#include <cvd/cpu_hacks.h>
#include <cvd/vector_image_ref.h>

#include <TooN/SymEigen.h>

#include <gvars3/instances.h>

#include <iostream>
#include <iterator>
#include <algorithm>
#include <functional>
#include <numeric>
#include <set>
#include <cmath>
#include <tr1/tuple>

using namespace CVD;
using namespace GVars3;
using namespace std::tr1;
using namespace std;
using namespace TooN;


//Create a gaussian window which is normalized so that the elements sum up to 1.
Image<float> make_normalised_gaussian_window(int s, float sigmas=3)
{
	ImageRef size(s, s);
	Image<float> w(size, 0);
	Vector<2> centre = vec(size)/2;
	

	float r_edge = s/2.0;
	//exp(-r_edge^2/2*sigma^2) = exp(-sigmas*sigma*sigmas*sigma/2*sigma*sigma) = exp(-sigmas^2/2)
	//therefore
	//r_edge^2/sigma^2 = sigmas^2
	//so sigma = r_edge/sigmas
	
	float sigma = r_edge/sigmas;
	ImageRef p(0,0);

	do{
		Vector<2> rr = vec(p) - centre;
		float r2 = rr*rr;
		w[p] = exp(-r2 / (2*sigma*sigma));
	} while(p.next(size));

	//Normalize the window
	
	float mul = accumulate(w.begin(), w.end(), 0.0);

	mul = 1/mul;

	transform(w.begin(), w.end(), w.begin(), bind1st(multiplies<float>(),mul));

	return w;
}




//Measure the similarity between a and b using weighted NCC. It is assumed that b
//has already been weighted and normalized, so only weighting and normalization on a
//needs to be computed. Furthermore, it is assumed that mask sums to 1, so computing the
//mean value of weighted a is easy.
//
//Given inputs with correct assumptions, this will produce the same results as ncc_compare_weighted
float ncc_compare_weighted_norm_b(const SubImage<float>& a, const SubImage<float>& b, const float Nb, const SubImage<float>& w)
{
	assert(a.size() == b.size());
	assert(a.size() == w.size());

	//w ***must*** sum to 1.
	float a_mean=0;
	//Compute the weighted mean and of weighted a.
	for(int y=0; y < a.size().y; y++)
		for(int x=0; x < a.size().x; x++)
			a_mean += a[y][x] * w[y][x];

	float ncc_top=0;
	float ncc_bot_a=0;
	for(int y=0; y < a.size().y; y++)
		for(int x=0; x < a.size().x; x++)
		{
			ncc_top   += (a[y][x] - a_mean)*b[y][x];
			ncc_bot_a += w[y][x] * (a[y][x] - a_mean) * (a[y][x] - a_mean);
		}

	return ncc_top / sqrt(ncc_bot_a * Nb);
}

float ncc_compare_weighted(const SubImage<float>& a, const SubImage<float>& b, const SubImage<float>& w)
{
	assert(a.size() == b.size());
	assert(a.size() ==w.size());

	//w ***must*** sum to 1.
	float a_mean=0, b_mean=0;
	//Compute the weighted mean of a and b
	for(int y=0; y < a.size().y; y++)
		for(int x=0; x < a.size().x; x++)
		{
			a_mean += a[y][x] * w[y][x];
			b_mean += b[y][x] * w[y][x];
		}

	float ncc_top=0;
	float ncc_bot_a=0;
	float ncc_bot_b=0;
	for(int y=0; y < a.size().y; y++)
		for(int x=0; x < a.size().x; x++)
		{
			ncc_top   += w[y][x] * (a[y][x] - a_mean) * (b[y][x] - b_mean);
			ncc_bot_a += w[y][x] * (a[y][x] - a_mean) * (a[y][x] - a_mean);
			ncc_bot_b += w[y][x] * (b[y][x] - b_mean) * (b[y][x] - b_mean);
		}

	return ncc_top / sqrt(ncc_bot_a * ncc_bot_b);
}

struct SortNCCLargestFirst
{
	bool operator()(const pair<float, Vector<2> >& a, const pair<float, Vector<2> >& b) const
	{
		return a.first > b.first;
	}
};

//Scan tmplte over image, and find the large local maxima, using weighted NCC.
//At most n are kept, and only if they exceed a threshold. mask must sum to 1.
vector<pair<float, Vector<2> > > ncc_scan(const SubImage<float>&image, const SubImage<float>& tmplte, const Image<float>& mask, unsigned int keep_n)
{
	assert(tmplte.size() == mask.size());
	ImageRef template_size = tmplte.size();
	ImageRef window_size = image.size() - template_size;
	
	//Normalize the template
	Image<float> normal_template(template_size);
	float norm_factor=0;
	{
		normal_template.copy_from(tmplte);
		//mask must sum to 1
		float mean=0;

		for(int y=0; y < tmplte.size().y; y++)
			for(int x=0; x < tmplte.size().x; x++)
				mean += tmplte[y][x] * mask[y][x];

		for(int y=0; y < tmplte.size().y; y++)
			for(int x=0; x < tmplte.size().x; x++)
			{
				normal_template[y][x] = (tmplte[y][x] - mean) * mask[y][x];
				norm_factor += mask[y][x] * (tmplte[y][x] - mean) * (tmplte[y][x] - mean);
			}
	}
	
	//Find the NCC all over the window, and store them if they have
	//a positive correlation
	Image<float> nccs(window_size, 0);

	for(int y=0; y < window_size.y; y++)
	{
		for(int x=0; x < window_size.x; x++)
		{
			ImageRef p(x, y);
			SubImage<float> si = image.sub_image(p, template_size);
			nccs[y][x] = ncc_compare_weighted_norm_b(si, normal_template, norm_factor, mask);
		}
	}
	
	//Keep the N best scores. There are not likely to be that many
	//maxima, so we can use a simple technique.
	vector<pair<float, Vector<2> > > maxima;

	for(int y=1; y < window_size.y-1; y++)
		for(int x=1; x < window_size.x-1; x++)
		{
			//Checkif it is a local maximum
			float ncc = nccs[y][x];

			if( ncc > nccs[y-1][x-1] &&
			    ncc > nccs[y-1][x  ] &&
			    ncc > nccs[y-1][x+1] &&
			    ncc > nccs[y  ][x-1] &&
			    ncc > nccs[y  ][x+1] &&
			    ncc > nccs[y+1][x-1] &&
			    ncc > nccs[y+1][x  ] &&
			    ncc > nccs[y+1][x+1] &&
				ncc > 0)
			{
				ImageRef p(x, y);
				pair<Vector<2>, double> extremum = interpolate_extremum_value(nccs, p);

				maxima.push_back(make_pair(extremum.second, extremum.first));
			}
		}
	
	//Sort the maxima (largest first), and keep the N largest
	sort(maxima.begin(), maxima.end(), SortNCCLargestFirst());
	if(maxima.size() > keep_n)
		maxima.resize(keep_n);

	return maxima;
}

int main(int argc, char** argv)
{
try
{	
	enableFPE();
	GUI.LoadFile("mosaic.cfg");
	int lastarg = GUI.parseArguments(argc, argv);

	if(lastarg != argc)
	{
		cerr << "Incorrect arguments!\n";
		return 1;
	}
	
	int tlen=GV3::get<int>("fine.mask.size", 0, -1);
	int slen=GV3::get<int>("fine.mask.search_radius", 0, -1);

	ImageRef template_size(tlen, tlen);
	ImageRef search_radius(slen, slen);
	ImageRef window_size = template_size + 2*search_radius;

	unsigned int keep_n = GV3::get<unsigned int>("fine.keep_n", 0, -1);
	
	//In order to overlay the input with the source, shift input by this value
	//So, pixel 0 of the input corresponds to pixel -offset of the image

	Vector<2> initial_offset = GV3::get<Vector<2> >("initial_offset", Zeros, -1);

	//This mask is used to smooth out the NCC, to prevent edge
	//effects.
	Image<float> mask = make_normalised_gaussian_window(tlen, GV3::get<float>("fine.mask.sigmas", 0, -1));

	//Load the image to be registered to
	Image<float> source(img_load(GV3::get<string>("source", "", -1)));
	
	//Load the image to be registered
	Image<float> image = img_load(GV3::get<string>("image", "", -1));


	//Use NCC to find an offset relative to ir_rounded(initial_offset)

	//For all possible subwindows in source which can be scanned
	for(int y=0; y < source.size().y - window_size.y; y++)
	{
		for(int x=0; x<source.size().x - window_size.x; x++)
		{
			//The window is the size of the template, plus a bit to allow for local
			//search
			ImageRef window_pos(x, y);
			ImageRef template_source_pos = window_pos + search_radius;
			ImageRef template_pos = template_source_pos- ir_rounded(initial_offset);

			//Snip a template out of the middle...
			
			if(image.in_image(template_pos) && image.in_image(template_pos+template_size))
			{
				SubImage<float> tmplte = image.sub_image(template_pos, template_size);

				vector<pair<float, Vector<2> > > best_nccs = ncc_scan(source.sub_image(window_pos, window_size), tmplte, mask, keep_n);
				
				//Aggregate together all the NCC scores and locations
				for(unsigned int i=0; i < best_nccs.size(); i++)
				{
					//The offset is caused by scanning the template over the window.
					//(0,0) for the window is at search_radius.
					cout << best_nccs[i].first << " " << best_nccs[i].second - vec(search_radius) + vec(ir_rounded(initial_offset)) << endl;
				}
			}

		}
	}
}
catch(Exceptions::All a)
{
	cerr << "Error: " << a.what << endl;
}
}
