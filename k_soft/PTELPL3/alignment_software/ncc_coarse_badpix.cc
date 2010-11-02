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
// and image patch b. 
float ncc_compare(const SubImage<float>& a, float a_mean, float std_a, const SubImage<float>& b, float b_mean, float std_b)
{
	assert(a.size() == b.size());

	float ncc=0;
	for(int y=0; y < a.size().y; y++)
		for(int x=0; x < a.size().x; x++)
			ncc += (a[y][x] - a_mean)*(b[y][x] - b_mean);

	//float std_a = sqrt(a_mean_sq - a_mean*a_mean);
	return ncc/(std_a*std_b);
}


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
vector<pair<float, ImageRef> > ncc_scan(const SubImage<float>&image, const SubImage<float>& tmplte, unsigned int keep_n, const Image<vector<pair<ImageRef, ImageRef> > >& bad_pixels)
{
	ImageRef template_size = tmplte.size();
	ImageRef window_size = image.size() - template_size;
	
	//Compute numbers required to normalize the template
	Image<float> tmplte_sq(tmplte.size());
	float template_sum=0, template_sum_sq=0;
	for(int r=0; r < tmplte.size().y; r++)
		for(int c=0; c < tmplte.size().x; c++)
		{
			template_sum+=tmplte[r][c];
			tmplte_sq[r][c] = tmplte[r][c] * tmplte[r][c];
			template_sum_sq+= tmplte_sq[r][c];
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
			
			//Extract a sub-image
			SubImage<float> si = image.sub_image(p, template_size);
			
			//Compute the normalization constants for the patch
			float si_sum = area_sum(integral, p, template_size);
			float si_sum_sq= area_sum(integral_sq, p, template_size);

			float template_sum_missing = template_sum;
			float template_sum_sq_missing = template_sum_sq;
			
			//Subtract off the missing pixels.
			int area = template_size.area() - bad_pixels[y][x].size();
			for(unsigned int m=0; m < bad_pixels[y][x].size(); m++)
			{
				assert(&image[bad_pixels[y][x][m].first] == &si[bad_pixels[y][x][m].second]);

				si_sum   -=image[bad_pixels[y][x][m].first];
				si_sum_sq-=squared_image[bad_pixels[y][x][m].first];
				
				template_sum_missing -= tmplte[bad_pixels[y][x][m].second];
				template_sum_sq_missing -= tmplte_sq[bad_pixels[y][x][m].second];
			}

			float si_mean = si_sum/area;
			float si_mean_sq = si_sum_sq / area;

			float template_mean = template_sum_missing / area;
			float template_mean_sq = template_sum_sq_missing / area;


			float si_std = sqrt(si_mean_sq - si_mean*si_mean);
			float template_std = sqrt(template_mean_sq - template_mean*template_mean);


			
			float ncc = ncc_compare(si, si_mean, si_std, tmplte, template_mean, template_std);

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

//Find the pixels marked as bad within a region, and store both the image and region
//relative offsets.
Image<vector<pair<ImageRef, ImageRef> > > missing_pixels(const Image<bool>& im, ImageRef w)
{
	Image<vector<pair<ImageRef, ImageRef> > > out(im.size() - w);
	for(int r=0; r < out.size().y; r++)
		for(int c=0; c < out.size().x; c++)
			for(int y=0; y < w.y; y++)
				for(int x=0; x < w.x; x++)
					if(!im[r+y][c+x])
						out[r][c].push_back(make_pair(ImageRef(x+c, r+y), ImageRef(x, y)));
	return out;
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

		//Load the badpixels image, if it exists.
		Image<vector<pair<ImageRef, ImageRef> > > bad_pixels(source.size());
		string bad_name = GV3::get<string>("badpixels", "", 1);

		if(bad_name != "")
			bad_pixels = missing_pixels(img_load(bad_name), tsize);
	

		unsigned int keep_n = GV3::get<unsigned int>("coarse.keep_n", 0, -1);
		float rq_radius = GV3::get<float>("coarse.case_deletion.radius", 0, -1);
		float rq_fraction = GV3::get<float>("coarse.case_deletion.fraction", 0, -1);

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
				vector<pair<float, ImageRef> > best_nccs = ncc_scan(source, tmplte, keep_n, bad_pixels);
				
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
	catch(Exceptions::All a)
	{
		cerr << "Error: " << a.what << endl;
	}
}
