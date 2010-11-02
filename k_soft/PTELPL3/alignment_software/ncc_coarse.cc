////////////////////////////////////////////////////////////////////////////////
//
// Astronomical image alignment software. Copyright Edward Rosten 2009.
//

#include <cvd/image_io.h>
#include <cvd/cpu_hacks.h>
#include <cvd/integral_image.h>
#include <cvd/draw.h>
#include <cvd/image_convert.h>
#include <cvd/convolution.h>
#include <cvd/vector_image_ref.h>

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

// This computes the NCC (normalized cross correlation) between image patch a
// and image patch b. The normalization factors are provided for a, and it 
// is assumed that b is normalized (ie mean 0, standard. dev 1).
float ncc_compare_norm_b(const SubImage<float>& a, float a_mean, float std_a, const SubImage<float>& b)
{
	assert(a.size() == b.size());

	float ncc=0;
	for(int y=0; y < a.size().y; y++)
		for(int x=0; x < a.size().x; x++)
			ncc += (a[y][x] - a_mean)*b[y][x];

	//float std_a = sqrt(a_mean_sq - a_mean*a_mean);
	return ncc/std_a;
}

// Computes the sum of all pixels over an area in an image. This is performed 
// using an integral image so it's efficient.
float area_sum(const SubImage<float>& integ, ImageRef pos, ImageRef size)
{
	if(pos == ImageRef(0,0))
		return integ[size - ImageRef(1,1)];
	else
	{
		int top = pos.y-1;
		int left = pos.x-1;
		
		int bottom = top + size.x;
		int right = left + size.y;

		if(pos.y == 0)//Ignore top row
			return integ[bottom][right] - integ[bottom][left];
		else if(pos.x == 0) //Ignore left column
			return integ[bottom][right] - integ[top][right];
		else
			return integ[top][left] + integ[bottom][right] - integ[top][right] - integ[bottom][left];
	}
}


// This function scans tmplte over image using NCC and returns the keep_n
// best NCC scores along with their positions. 
vector<pair<float, ImageRef> > ncc_scan(const SubImage<float>&image, const SubImage<float>& tmplte, unsigned int keep_n)
{
	ImageRef template_size = tmplte.size();
	ImageRef window_size = image.size() - template_size;
	
	//Compute a normalized template and store the result in normal_template
	Image<float> normal_template;
	{
		normal_template.copy_from(tmplte);

		float mean=0, mean_sq=0;
		for(Image<float>::iterator i=normal_template.begin(); i != normal_template.end(); i++)
		{
			mean+=*i;
			mean_sq+=*i**i;
		}
		mean/=normal_template.size().area();
		mean_sq/=normal_template.size().area();
		float inv_std = 1/sqrt(mean_sq - mean*mean);
		
		//Ensure the normalized template is zero mean, unit variance.
		for(Image<float>::iterator i=normal_template.begin(); i != normal_template.end(); i++)
			*i = (*i-mean)*inv_std;
	}

	//Create integral images for fast computation of the 
	//normalization factors for all subimages.
	Image<float> integral = integral_image(image);

	Image<float> squared_image;
	squared_image.copy_from(image);
	for(Image<float>::iterator i=squared_image.begin(); i != squared_image.end(); i++)
		*i*=*i;
	
	Image<float> integral_sq = integral_image(squared_image);

	//Keep the N best ncc scores, using a min-heap. This allows us to always
	//remove the smallest element, keeping the largest ones.

	//C++ heap functions use std::less to create a max-heap, growing from the
	//beginning of the array. This allows for convenient sorting. pop_heap
	//gets the largest element, and the heap is shrunk by 1. This element
	//is then put on to the end of the array, ie the spot freed up by shrinking
	//the heap by 1. Repeating this procedure will sort the heap, pulling out
	//the largest elements and placing them at the end. The resulting array will
	//then be sorted by std::less.

	//Therefore we need to use std::greater to create a min-heap

	//The first element in the array will be the smallest value
	typedef greater<pair<float, ImageRef> > minheap_compare;

	vector<pair<float, ImageRef> > ncc_heap;

	for(int y=0; y < window_size.y; y++)
	{
		for(int x=0; x < window_size.x; x++)
		{
			ImageRef p(x, y);

			SubImage<float> si = image.sub_image(p, template_size);

			float si_mean = area_sum(integral, p, template_size)/template_size.area();
			float si_mean_sq= area_sum(integral_sq, p, template_size)/template_size.area();
			float si_std = sqrt(si_mean_sq - si_mean*si_mean);
			
			float ncc = ncc_compare_norm_b(si, si_mean, si_std, normal_template);

			if(ncc_heap.size() < keep_n || ncc_heap.size() == 0 || ncc > ncc_heap[0].first)
			{
				ncc_heap.push_back(make_pair(ncc, p));
				push_heap(ncc_heap.begin(), ncc_heap.end(), minheap_compare());
			}

			if(ncc_heap.size() > keep_n)
			{
				pop_heap(ncc_heap.begin(), ncc_heap.end(), minheap_compare());
				ncc_heap.pop_back();
			}
		}
	}
	return ncc_heap;
}

int main(int argc, char** argv)
{
	try
	{	
		enableFPE();
		GUI.LoadFile("mosaic.cfg");
		int lastarg = GUI.parseArguments(argc, argv);
		
		//Load all the unregistered images
		vector<string> image_names;
		for(int i=lastarg; i < argc; i++)
			image_names.push_back(argv[i]);

		//Load the image to which they will all be registered
		string source_name = GV3::get<string>("source", "", -1);
		Image<float> source(img_load(source_name));

		int tlen=GV3::get<int>("coarse.chunk.size", 0, -1);
		int num_chunks=GV3::get<int>("coarse.chunk.num_per_side", 0, -1);
		int border_len=GV3::get<int>("coarse.chunk.border", 0, -1);

		ImageRef tsize(tlen, tlen), border(border_len, border_len);

		unsigned int keep_n = GV3::get<unsigned int>("coarse.keep_n", 0, -1);
		float rq_radius = GV3::get<float>("coarse.case_deletion.radius", 0, -1);
		float rq_fraction = GV3::get<float>("coarse.case_deletion.fraction", 0, -1);

		string method = GV3::get<string>("coarse.method", "", -1);
		float histogram_radius = GV3::get<float>("coarse.histogram.radius", 0, -1);

		Image<Rgb<byte> > out(source.size()*3, Rgb<byte>(0,0,255));
		ImageRef out_off=source.size();

		//For the image, snip out a chunk and scan it all over the source image
		//to find a registration for the chunk. Repeat for many chunks.
		
		string image_name = GV3::get<string>("image", "", -1);
		Image<float> image = img_load(image_name);

		vector<ImageRef> ncc_offsets;

		ImageRef best_pos, best_off;
		
		for(int ny=0; ny<=num_chunks; ny++)
			for(int nx=0; nx<=num_chunks; nx++)
			{
				ImageRef last_pos = image.size() - tsize - 2*border;

				ImageRef pos; 

				if(num_chunks == 1)
					pos = (image.size() - tsize)/2;
				else
					pos = border + ImageRef(last_pos.x*nx/(num_chunks-1), last_pos.y*ny/(num_chunks-1));

				//Snip a template out of the middle...
				SubImage<float> tmplte = image.sub_image(pos, tsize);
				
				//Scan the template all over the source image:
				vector<pair<float, ImageRef> > best_nccs = ncc_scan(source, tmplte, keep_n);
				
				//Aggregate together all the NCC scores and locations
				for(unsigned int i=0; i < best_nccs.size(); i++)
				{
					//float ncc = best_nccs[i].first;
					ImageRef off = best_nccs[i].second - pos;
					ncc_offsets.push_back(off);
				}
			}
		


		//At this point we have scanned all chunks of image over the source
		///image.

		if(method == "histogram")
		{
			ImageRef h_size = source.size() * 3;
			ImageRef h_off = h_size / 2;
			
			//Build up a histogram of offsets
			Image<float> histogram(h_size, 0);
			for(unsigned int i=0; i < ncc_offsets.size(); i++)
				histogram[h_off + ncc_offsets[i]]++;
			
			//Blur the histogram a bit
			convolveGaussian(histogram, histogram, histogram_radius);

			//Find the high point.

			ImageRef offset = histogram.pos(max_element(histogram.begin(), histogram.end())) - h_off;
			cout << vec(offset) << endl;
		}
		else
		{
			Vector<2> mean;
			//Perform some case deletion robust estimation to find the real offset
			for(int itnum=0; ; itnum++)
			{
				//Compute mean offset position
				ImageRef mean_sum(0,0);

				for(unsigned int i=0; i < ncc_offsets.size(); i++)
					mean_sum += ncc_offsets[i];

				mean = vec(mean_sum) / ncc_offsets.size();

				//Keep the best N% of the points: first record the error for each point.
				//Also record the covariance for the termination condition.
				vector<pair<float, int> > radius_and_index;
				Matrix<2> cov = Zeros;
				
				for(unsigned int i=0; i < ncc_offsets.size(); i++)
				{
					//Compute the covariance
					Vector<2> d = vec(ncc_offsets[i]) - mean;
					cov += d.as_col() * d.as_row();

					//Record the errors.
					radius_and_index.push_back(make_pair((float)(d*d), i));
				}

				//Terminate when the largest spread in the covariance
				//drops below some threshold
				cov /= ncc_offsets.size();
				SymEigen<2> ec(cov);
				float sig = sqrt(ec.get_evalues()[1]);
				if(sig  < rq_radius)
					break;
			
				
				//Sort the points so that the best ones (smallest radius) are at the
				//beginning of the array
				sort(radius_and_index.begin(), radius_and_index.end());

				//Discard the indices of the worst points.
				int num_to_keep = (int)floor(rq_fraction * radius_and_index.size());
				radius_and_index.resize(num_to_keep);
				
				//Keep only the best points.
				vector<ImageRef> ncc_offsets_shrunk;
				for(unsigned int i=0; i < radius_and_index.size(); i++)
					ncc_offsets_shrunk.push_back(ncc_offsets[radius_and_index[i].second]);

				ncc_offsets = ncc_offsets_shrunk;
			}
			
			cout << mean << endl;
		}
	}
	catch(Exceptions::All a)
	{
		cerr << "Error: " << a.what << endl;
	}
}
