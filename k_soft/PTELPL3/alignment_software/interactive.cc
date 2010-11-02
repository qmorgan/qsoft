////////////////////////////////////////////////////////////////////////////////
//
// Astronomical image alignment software. Copyright Edward Rosten 2009.
//

#include <cvd/image_io.h>
#include <cvd/image_convert.h>
#include <cvd/connected_components.h>
#include <cvd/glwindow.h>
#include <cvd/draw.h>
#include <cvd/gl_helpers.h>

#include <gvars3/instances.h>
#include <gvars3/GUI_readline.h>
#include <gvars3/GStringUtil.h>

#include <iostream>
#include <map>
#include <string>
#include <iterator>
#include <algorithm>
#include <math.h>

using namespace CVD;
using namespace GVars3;
using namespace std;


map<string, string> watch;

void quit(void*, string comm, string d)
{
	GV3::set_var("keep_running", "0");
}
void watch_var(void*, string comm, string d)
{
	vector<string> vs = ChopAndUnquoteString(d);

	if(vs.size() != 1)
	{
		cerr << "Error: " << comm << " takes 1 argument: " << comm << " gvar\n";
		return;
	}

	watch[vs[0]] = GV3::get_var(vs[0]);
}

bool watch_update()
{
	bool changes=0;

	for(map<string, string>::iterator i=watch.begin(); i != watch.end(); i++)
	{
		string s = GV3::get_var(i->first);

		if(s != watch[i->first])
		{
			changes=1;
			watch[i->first] = s;
		}
	}
	
	return changes;
}



float lim_pow(float f, float p)
{
	if(f < 0)
		return 0;
	else
	{
		f = pow(f, p);
		if(f > 1)
			return 1;
		else 
			return f;
	}
}

int main(int argc, char** argv)
{
	try
	{	
		GUI.RegisterCommand("watch", watch_var);
		GUI.RegisterCommand("quit", quit);
		GUI.RegisterCommand("exit", quit);

		int lastarg = GUI.parseArguments(argc, argv);

		vector<string> image_names;
		for(int i=lastarg; i < argc; i++)
			image_names.push_back(argv[i]);

		vector<Image<float> > images;
		vector<vector<float> > sorted_pixels;
		for(unsigned int i=0;i <  image_names.size(); i++)
		{
			images.push_back(img_load(image_names[i]));

			vector<float> srt(images.back().begin(), images.back().end());
			sort(srt.begin(), srt.end());
			sorted_pixels.push_back(srt);
		}

		sorted_pixels.resize(images.size());

		Image<float> i = images.at(0);

		float minval = *min_element(i.begin(), i.end());
		float maxval = *max_element(i.begin(), i.end());

		GV3::get<float>("max",0,1) = maxval*2;
		GV3::get<float>("nimages",0,1) = images.size();
		gvar3<float> centre("centre", (maxval+minval)/2, 1);
		gvar3<float> range("range", maxval-minval, 1);
		gvar3<int> imagenum("imagenum", 0, 1);
		gvar3<float> percentile("percentile", 0, 1);
		gvar3<float> radius("radius", 0, 1);
		gvar3<int> min_blob_size("min_blob_size", 0, 1);

		GUI.LoadFile("interactive.cfg");

		ImageRef imsize = i.size();

		GLWindow d(i.size());
		float zoom=1;

		Image<byte> s(i.size());
		Image<Rgb<byte> > out(i.size());
		Image<byte> deadpix(i.size());

		readline_in_current_thread line("> ");

		gvar3<bool> keep_running("keep_running", 1, 1), redraw("redraw", 1, 1);
		gvar3<int> invert("invert", 0, 1);
		bool had_update=0;

		vector<float>* sp = &sorted_pixels[0];
		
		while(*keep_running)
		{
			GUI_Widgets.process_in_crnt_thread();
			line.poll();

			if(*imagenum < (int)images.size() && *imagenum >= 0)
			{
				i = images[*imagenum];
				sp = &sorted_pixels[*imagenum];
				GV3::get<string>("image_name","",1) = argv[*imagenum+lastarg];
			}
			if(d.has_events())
			{
				vector<GLWindow::Event> e;
				d.get_events(e);

				for(unsigned int i=0; i < e.size(); i++)
				{
					if(e[i].type == GLWindow::Event::RESIZE)
					{
						ImageRef newsize = e[i].size;
						ImageRef size=imsize;

						float old_r = (float)imsize.x / imsize.y;
						float new_r = (float)newsize.x / newsize.y;


						glViewport(0, 0, newsize.x, newsize.y);
						glMatrixMode(GL_PROJECTION);
						glLoadIdentity();

						if(new_r > old_r) //Then use the y axis
							zoom = newsize.y / (float)size.y; 
						else
							zoom = newsize.x / (float)size.x;

						glOrtho(-.5/zoom, (newsize.x-1.5)/zoom, (newsize.y-1.5)/zoom, -.5/zoom, -1 , 1);

						glPixelZoom(zoom,-zoom);
						glRasterPos2f(0, 0);
						*redraw=1;
					}
					else if(e[i].type == GLWindow::Event::EVENT)
					{
						if(e[i].which == GLWindow::EVENT_CLOSE)
							*keep_running=0;
						else if(e[i].which == GLWindow::EVENT_EXPOSE)
							*redraw = 1;
					}
				}

			}
			
			if(watch_update())
				had_update=1;

			if(*redraw || (had_update && !watch_update()))
			{
				float ctr = *centre, rng = *range, gm = GV3::get<float>("gamma", 1, 1);
			
				for(int r=0; r < i.size().y; r++)
					for(int c=0; c < i.size().y; c++)
						Pixel::DefaultConversion<float, byte>::type::convert(lim_pow((i[r][c] - ctr)/rng, gm), s[r][c]);
				if(*invert)
					for(int r=0; r < i.size().y; r++)
						for(int c=0; c < i.size().y; c++)
							s[r][c] = 255 - s[r][c];

				float threshold = (*sp)[(int)floor(sp->size() * *percentile)];
				
				out = convert_image(s);

				vector<ImageRef> d_deadpix;
				deadpix.fill(0);
				vector<ImageRef> disc=getDisc(*radius);

				for(int r=0; r < i.size().y; r++)
					for(int c=0; c < i.size().y; c++)
						if(i[r][c] < threshold)
						{
							ImageRef p(c,r);
							deadpix[p]=1;
							for(unsigned int j=0; j < disc.size(); j++)
								d_deadpix.push_back(p + disc[j]);
						}

				vector<vector<ImageRef> > bits, ibits;
				connected_components(d_deadpix, bits);

				ibits.resize(bits.size());
				//Intersect
				for(unsigned int i=0; i < bits.size(); i++)
					for(unsigned int j=0; j < bits[i].size(); j++)
						if(deadpix.in_image(bits[i][j]) && deadpix[bits[i][j]])
							ibits[i].push_back(bits[i][j]);

				for(unsigned int i=0; i < ibits.size(); i++)
					for(unsigned int j=0; j < ibits[i].size(); j++)
						if((int)ibits[i].size() < *min_blob_size)
							out[ibits[i][j]].red=255;
						else
							out[ibits[i][j]].green=255;
						
						

				glDrawPixels(out);
				glFlush();
				d.swap_buffers();
				*redraw=0;
				had_update=0;
			}
			usleep(100000);

		}
	}
	catch(Exceptions::All a)
	{
		cerr << "Error: " << a.what << endl;
	}
}
