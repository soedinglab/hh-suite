//**********  pngwriter.cc   **********************************************
//  Author:                    Paul Blackburn
//
//  Email:                     individual61@users.sourceforge.net
//
//  Version:                   0.3.6   (18 / IX / 2003)
//
//  Description:               Object that allows plotting a 48 bit
//                             PNG image pixel by pixel, which can 
//                             then be opened with a graphics program.
//  
//  License:                   GNU General Public License
//                             © 2002, 2003 Paul Blackburn
//                             
//  Website: Main:             http://pngwriter.sourceforge.net/
//           Sourceforge.net:  http://sourceforge.net/projects/pngwriter/
//           Freshmeat.net:    http://freshmeat.net/projects/pngwriter/
//           
//  Documentation:             The header file (pngwriter.h) is commented, but for a
//                             quick reference document, and support,
//                             take a look at the website.
//
//*************************************************************************

#include "pngwriter.h"

//Constructor for int colour levels, char * filename
//////////////////////////////////////////////////////////////////////////
pngwriter::pngwriter(int x, int y, int backgroundcolour, char * filename)
{
   width_ = x;
   height_ = y;
   backgroundcolour_ = backgroundcolour;
   compressionlevel_ = -2;
   filegamma_ = 0.5;
   
   strcpy(textauthor_, "PNGwriter Author: Paul Blackburn");
   strcpy(textdescription_, "http://pngwriter.sourceforge.net/");
   strcpy(textsoftware_, "PNGwriter: a C++ class to make plotting to PNG files easy.");
   strcpy(texttitle_, filename);  

   if((width_<0)||(height_<0))
     {
	cerr << "PNGwriter::pngwriter - ERROR **: Constructor called with negative height or width." << endl;
     }
   
   
   if(backgroundcolour_ >65535)
     {
	cerr << "PNGwriter::pngwriter - WARNING **: Constructor called with background colour greater than 65535. Setting to 65535."<<endl;
     }
   
   if(backgroundcolour_ <0)
     {
	cerr << "PNGwriter::pngwriter - WARNING **: Constructor called with background colour lower than 0. Setting to 0."<<endl;
     }
   
   havesinecosinetables_ = 0; //Flag to avoid calculating sines and cosines over and over again.
   int kkkk;
   
   char * temp;
   
   bit_depth_ = 16; //Default bit depth for new images
   colortype_=2;
   screengamma_ = 2.2;
   
   temp = filename;
   
   strcpy(filename_, filename);
   
   graph_ = (png_bytepp)malloc(height_ * sizeof(png_bytep));   
   
   for (kkkk = 0; kkkk < height_; kkkk++)
     {
        graph_[kkkk] = (png_bytep)malloc(6*width_ * sizeof(png_byte));
	if(graph_[kkkk] == NULL)
	  {   
	     cerr << "PNGwriter::pngwriter - ERROR **:  Not able to allocate memory for image." << endl;
	  }
     }
  
   if(graph_ == NULL)
     {   
	cerr << "PNGwriter::pngwriter - ERROR **:  Not able to allocate memory for image." << endl;
     }
   
   for(int hhh = 0; hhh<width_;hhh++)
     {   
	for(int vhhh = 0; vhhh<height_;vhhh++)
	  {
	     
	     graph_[vhhh][6*hhh] = (char) floor(((double)backgroundcolour_)/256);
	     graph_[vhhh][6*hhh+1] = (char)(backgroundcolour_%256);
	     graph_[vhhh][6*hhh+2] = (char) floor(((double)backgroundcolour_)/256);
	     graph_[vhhh][6*hhh+3] = (char)(backgroundcolour_%256);
	     graph_[vhhh][6*hhh+4] = (char) floor(((double)backgroundcolour_)/256);
	     graph_[vhhh][6*hhh+5] = (char)(backgroundcolour_%256);
	  }
     }
   
};

//Constructor for double levels, char * filename
/////////////////////////////////////////////////////////////////////////
pngwriter::pngwriter(int x, int y, double backgroundcolour, char * filename)
{
   width_ = x;
   height_ = y;
   compressionlevel_ = -2;
   filegamma_ = 0.6;
   backgroundcolour_ = int(backgroundcolour*65535);
   
   strcpy(textauthor_, "PNGwriter Author: Paul Blackburn");
   strcpy(textdescription_, "http://pngwriter.sourceforge.net/");
   strcpy(textsoftware_, "PNGwriter: a C++ class to make plotting to PNG files easy.");
   strcpy(texttitle_, filename);   

   if((width_<0)||(height_<0))
     {
	cerr << "PNGwriter::pngwriter - ERROR **: Constructor called with negative height or width." << endl;
     }
   
   
   if(backgroundcolour_ >65535)
     {
	cerr << "PNGwriter::pngwriter - WARNING **: Constructor called with background colour greater than 1.0. Setting to 1.0."<<endl;
     }
   
   if(backgroundcolour_ <0)
     {
	cerr << "PNGwriter::pngwriter - WARNING **: Constructor called with background colour lower than 0.0. Setting to 0.0."<<endl;
     }
   
   havesinecosinetables_ = 0; //Flag to avoid calculating sines and cosines over and over again.
   int kkkk;
   
   char * temp;
   
   bit_depth_ = 16; //Default bit depth for new images
   colortype_=2;
   screengamma_ = 2.2;
   
   temp = filename;
   
   strcpy(filename_, filename);
   
   graph_ = (png_bytepp)malloc(height_ * sizeof(png_bytep));   
   
   for (kkkk = 0; kkkk < height_; kkkk++)
     {
        graph_[kkkk] = (png_bytep)malloc(6*width_ * sizeof(png_byte));
	if(graph_[kkkk] == NULL)
	  {   
	     cerr << "PNGwriter::pngwriter - ERROR **:  Not able to allocate memory for image." << endl;
	  }
     }
  
   if(graph_ == NULL)
     {   
	cerr << "PNGwriter::pngwriter - ERROR **:  Not able to allocate memory for image." << endl;
     }
   
   for(int hhh = 0; hhh<width_;hhh++)
     {   
	for(int vhhh = 0; vhhh<height_;vhhh++)
	  {
	     
	     graph_[vhhh][6*hhh] = (char) floor(((double)backgroundcolour_)/256);
	     graph_[vhhh][6*hhh+1] = (char)(backgroundcolour_%256);
	     graph_[vhhh][6*hhh+2] = (char) floor(((double)backgroundcolour_)/256);
	     graph_[vhhh][6*hhh+3] = (char)(backgroundcolour_%256);
	     graph_[vhhh][6*hhh+4] = (char) floor(((double)backgroundcolour_)/256);
	     graph_[vhhh][6*hhh+5] = (char)(backgroundcolour_%256);
	  }
     }
   
};

//Destructor
///////////////////////////////////////
pngwriter::~pngwriter()
{
   int jjj;
   for (jjj = 0; jjj < height_; jjj++) free(graph_[jjj]);            
   free(graph_);
   
};

//Constructor for int levels, const char * filename
//////////////////////////////////////////////////////////////
pngwriter::pngwriter(int x, int y, int backgroundcolour, const char * filename)
{
   width_ = x;
   height_ = y;
   backgroundcolour_ = backgroundcolour;
   compressionlevel_ = -2;
   filegamma_ = 0.6;
   
   strcpy(textauthor_, "PNGwriter Author: Paul Blackburn");
   strcpy(textdescription_, "http://pngwriter.sourceforge.net/");
   strcpy(textsoftware_, "PNGwriter: a C++ class to make plotting to PNG files easy.");
   strcpy(texttitle_, filename);   

   if((width_<0)||(height_<0))
     {
	cerr << "PNGwriter::pngwriter - ERROR **: Constructor called with negative height or width." << endl;
     }
   
   
   if(backgroundcolour_ >65535)
     {
	cerr << "PNGwriter::pngwriter - WARNING **: Constructor called with background colour greater than 65535. Setting to 65535."<<endl;
     }
   
   if(backgroundcolour_ <0)
     {
	cerr << "PNGwriter::pngwriter - WARNING **: Constructor called with background colour lower than 0. Setting to 0."<<endl;
     }
   
   havesinecosinetables_ = 0; //Flag to avoid calculating sines and cosines over and over again.
   int kkkk;
   
   char * temp;
   
   bit_depth_ = 16; //Default bit depth for new images
   colortype_=2;
   screengamma_ = 2.2;
   
   temp = (char *)filename;
   
   strcpy(filename_, filename);
   
   graph_ = (png_bytepp)malloc(height_ * sizeof(png_bytep));   
   
   for (kkkk = 0; kkkk < height_; kkkk++)
     {
        graph_[kkkk] = (png_bytep)malloc(6*width_ * sizeof(png_byte));
	if(graph_[kkkk] == NULL)
	  {   
	     cerr << "PNGwriter::pngwriter - ERROR **:  Not able to allocate memory for image." << endl;
	  }
     }
  
   if(graph_ == NULL)
     {   
	cerr << "PNGwriter::pngwriter - ERROR **:  Not able to allocate memory for image." << endl;
     }
   
   for(int hhh = 0; hhh<width_;hhh++)
     {   
	for(int vhhh = 0; vhhh<height_;vhhh++)
	  {
	     
	     graph_[vhhh][6*hhh] = (char) floor(((double)backgroundcolour_)/256);
	     graph_[vhhh][6*hhh+1] = (char)(backgroundcolour_%256);
	     graph_[vhhh][6*hhh+2] = (char) floor(((double)backgroundcolour_)/256);
	     graph_[vhhh][6*hhh+3] = (char)(backgroundcolour_%256);
	     graph_[vhhh][6*hhh+4] = (char) floor(((double)backgroundcolour_)/256);
	     graph_[vhhh][6*hhh+5] = (char)(backgroundcolour_%256);
	  }
     }
   
};

//Constructor for double int levels, const char * filename
/////////////////////////////////////////////////////////////////////////
pngwriter::pngwriter(int x, int y, double backgroundcolour, const char * filename)
{
   width_ = x;
   height_ = y;
   compressionlevel_ = -2;
   backgroundcolour_ = int(backgroundcolour*65535);
   filegamma_ = 0.6;
   
   strcpy(textauthor_, "PNGwriter Author: Paul Blackburn");
   strcpy(textdescription_, "http://pngwriter.sourceforge.net/");
   strcpy(textsoftware_, "PNGwriter: a C++ class to make plotting to PNG files easy.");
   strcpy(texttitle_, filename);   

   if((width_<0)||(height_<0))
     {
	cerr << "PNGwriter::pngwriter - ERROR **: Constructor called with negative height or width." << endl;
     }
   
   
   if(backgroundcolour_ >65535)
     {
	cerr << "PNGwriter::pngwriter - WARNING **: Constructor called with background colour greater than 1.0. Setting to 1.0."<<endl;
     }
   
   if(backgroundcolour_ <0)
     {
	cerr << "PNGwriter::pngwriter - WARNING **: Constructor called with background colour lower than 0.0. Setting to 0.0."<<endl;
     }
   
   havesinecosinetables_ = 0; //Flag to avoid calculating sines and cosines over and over again.
   int kkkk;
   
   char * temp;
   
   bit_depth_ = 16; //Default bit depth for new images
   colortype_=2;
   screengamma_ = 2.2;
   
   temp = (char *)filename;
   
   strcpy(filename_, filename);
   
   graph_ = (png_bytepp)malloc(height_ * sizeof(png_bytep));   
   
   for (kkkk = 0; kkkk < height_; kkkk++)
     {
        graph_[kkkk] = (png_bytep)malloc(6*width_ * sizeof(png_byte));
	if(graph_[kkkk] == NULL)
	  {   
	     cerr << "PNGwriter::pngwriter - ERROR **:  Not able to allocate memory for image." << endl;
	  }
     }
  
   if(graph_ == NULL)
     {   
	cerr << "PNGwriter::pngwriter - ERROR **:  Not able to allocate memory for image." << endl;
     }
   
   for(int hhh = 0; hhh<width_;hhh++)
     {   
	for(int vhhh = 0; vhhh<height_;vhhh++)
	  {
	     
	     graph_[vhhh][6*hhh] = (char) floor(((double)backgroundcolour_)/256);
	     graph_[vhhh][6*hhh+1] = (char)(backgroundcolour_%256);
	     graph_[vhhh][6*hhh+2] = (char) floor(((double)backgroundcolour_)/256);
	     graph_[vhhh][6*hhh+3] = (char)(backgroundcolour_%256);
	     graph_[vhhh][6*hhh+4] = (char) floor(((double)backgroundcolour_)/256);
	     graph_[vhhh][6*hhh+5] = (char)(backgroundcolour_%256);
	  }
     }
   
};


///////////////////////////////////////////////////////////////
void pngwriter::plot(int x, int y, int red, int green, int blue)
{
   
   if((bit_depth_ == 16))
     {	
	if( (height_-y >-1) && (height_-y <height_) && (6*(x-1) >-1) && (6*(x-1)+5<6*width_) )
	  {
	     graph_[height_-y][6*(x-1)] = (char) floor(((double)red)/256);
	     graph_[height_-y][6*(x-1)+1] = (char)(red%256);
	     graph_[height_-y][6*(x-1)+2] = (char) floor(((double)green)/256);
	     graph_[height_-y][6*(x-1)+3] = (char)(green%256);
	     graph_[height_-y][6*(x-1)+4] = (char) floor(((double)blue)/256);
	     graph_[height_-y][6*(x-1)+5] = (char)(blue%256);
	  };   	   
	
	/*
	 if(!( (height_-y >-1) && (height_-y <height_) && (6*(x-1) >-1) && (6*(x-1)+5<6*width_) ))
	 {
	 cerr << "PNGwriter::plot-- Plotting out of range! " << y << "   " << x << endl;
	 }
	 */   	
     }
   
   if((bit_depth_ == 8))
     {
	if( (height_-y >-1) && (height_-y <height_) && (3*(x-1) >-1) && (3*(x-1)+5<3*width_) )
	  {
	    
	     graph_[height_-y][3*(x-1)] = (char)(red%256);
	    
	     graph_[height_-y][3*(x-1)+1] = (char)(green%256);
	    
	     graph_[height_-y][3*(x-1)+2] = (char)(blue%256);
	  };   	   
	
	/*
	 if(!( (height_-y >-1) && (height_-y <height_) && (6*(x-1) >-1) && (6*(x-1)+5<6*width_) ))
	 {
	 cerr << "PNGwriter::plot-- Plotting out of range! " << y << "   " << x << endl;
	 }
	 */   
     }
};

void pngwriter::plot(int x, int y, double red, double green, double blue)
{
   plot(x,y,int(red*65535),int(green*65535),int(blue*65535));
};

///////////////////////////////////////////////////////////////
int pngwriter::read(int x, int y, int colour)
{
   int temp1;

   if((colour !=1)&&(colour !=2)&&(colour !=3))
     {
	cerr << "PNGwriter::read - WARNING **: Invalid argument: should be 1, 2 or 3, is " << colour << endl;
     }
   
   if(bit_depth_ == 16)
     {
	if(colour == 1)
	  {
	     temp1 = (graph_[(height_-y)][6*(x-1)])*256 + graph_[height_-y][6*(x-1)+1];
	     return temp1;
	  }
	
	if(colour == 2)
	  {
	     temp1 = (graph_[height_-y][6*(x-1)+2])*256 + graph_[height_-y][6*(x-1)+3];
	     return temp1;
	  }
	
	if(colour == 3)
	  {
	     temp1 = (graph_[height_-y][6*(x-1)+4])*256 + graph_[height_-y][6*(x-1)+5];
	     return temp1;
	  }
     }

   if(bit_depth_ == 8)
     {
	if(colour == 1)
	  {
	     temp1 = graph_[height_-y][3*(x-1)];
	     return temp1*256;
	  }
	
	if(colour == 2)
	  {
	     temp1 =  graph_[height_-y][3*(x-1)+1];
	     return temp1*256;
	  }
	
	if(colour == 3)
	  {
	     temp1 =  graph_[height_-y][3*(x-1)+2];
	     return temp1*256;
	  }
     }

   return 0; // should never get here
}


///////////////////////////////////////////////////////////////
int pngwriter::read(int x, int y)
{
   int temp1,temp2,temp3,temp4;
   
   if(bit_depth_ == 16)
     {
	temp1 = (graph_[(height_-y)][6*(x-1)])*256 + graph_[height_-y][6*(x-1)+1];
	temp2 = (graph_[height_-y][6*(x-1)+2])*256 + graph_[height_-y][6*(x-1)+3];
	temp3 = (graph_[height_-y][6*(x-1)+4])*256 + graph_[height_-y][6*(x-1)+5];
	temp4 =  int((temp1+temp2+temp3)/3.0);  
     }
   else if(bit_depth_ == 8)
     {
	temp1 = graph_[height_-y][3*(x-1)];
	temp2 =  graph_[height_-y][3*(x-1)+1];
	temp3 =  graph_[height_-y][3*(x-1)+2];
	temp4 =  int((temp1+temp2+temp3)/3.0);
     }
   else
     {
	cerr << "PNGwriter::read - WARNING **: Invalid bit depth! Returning 0 as average value." << endl;
	temp4 = 0;
     }

   return temp4;
   
}

/////////////////////////////////////////////////////
double  pngwriter::dread(int x, int y, int colour)
{
   return double(read(x,y,colour))/65535.0;   
}

double  pngwriter::dread(int x, int y)
{
   return double(read(x,y))/65535.0;   
}

///////////////////////////////////////////////////////   
void pngwriter::clear()
{
   int pen = 0;
   int pencil = 0;
   
   if(bit_depth_==16)
     {
	for(pen = 0; pen<width_;pen++)
	  {
	     for(pencil = 0; pencil<height_;pencil++)
	       {
		  
		  graph_[pencil][6*pen] = 0;
		  graph_[pencil][6*pen+1] = 0;
		  graph_[pencil][6*pen+2] = 0;
		  graph_[pencil][6*pen+3] = 0;
		  graph_[pencil][6*pen+4] = 0;
		  graph_[pencil][6*pen+5] = 0;	     
	       }
	  }
     }   
   
   if(bit_depth_==8)
     {
	for(pen = 0; pen<width_;pen++)
	  {
	     for(pencil = 0; pencil<height_;pencil++)
	       {
		  graph_[pencil][3*pen] = 0;
		  graph_[pencil][3*pen+1] = 0;
		  graph_[pencil][3*pen+2] = 0;
	       }
	  }
     }
   
};

/////////////////////////////////////////////////////
void pngwriter::pngwriter_rename(char * newname)
{	     
   strcpy(filename_,newname);
   strcpy(texttitle_,newname);
};


///////////////////////////////////////////////////////
void pngwriter::pngwriter_rename(const char * newname)
{	     
   strcpy(filename_,newname);
   strcpy(texttitle_,newname);   
};

///////////////////////////////////////////////////////
void pngwriter::settext(char * title, char * author, char * description, char * software)
{
  strcpy(texttitle_, title);
  strcpy(textauthor_, author);
  strcpy(textdescription_, description);
  strcpy(textsoftware_, software);
};

///////////////////////////////////////////////////////
void pngwriter::settext(const char * title, const char * author, const char * description, const char * software)
{
  strcpy(texttitle_, title);
  strcpy(textauthor_, author);
  strcpy(textdescription_, description);
  strcpy(textsoftware_, software);
};

///////////////////////////////////////////////////////
void pngwriter::close()
{
   FILE            *fp;
   png_structp     png_ptr;
   png_infop       info_ptr;
   
   fp = fopen(filename_, "wb");                            
   png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
   info_ptr = png_create_info_struct(png_ptr);             
   png_init_io(png_ptr, fp);
   if(compressionlevel_ != -2)
     {
	png_set_compression_level(png_ptr, compressionlevel_);                               
     }
   else
     {
	png_set_compression_level(png_ptr, PNGWRITER_DEFAULT_COMPRESSION);                               
     }
	
   png_set_IHDR(png_ptr, info_ptr, width_, height_,         
		bit_depth_, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

   if(filegamma_ < 1.0e-1)
     {
	filegamma_ = 0.7;
     }
   
   png_set_gAMA(png_ptr, info_ptr, filegamma_);

   time_t          gmt;
   png_time        mod_time;
   png_text        text_ptr[5];
   time(&gmt);
   png_convert_from_time_t(&mod_time, gmt);
   png_set_tIME(png_ptr, info_ptr, &mod_time);
   text_ptr[0].key = new(char[6]); strcpy(text_ptr[0].key,"Title");
   text_ptr[0].text = texttitle_;
   text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
   text_ptr[1].key = new(char[7]); strcpy(text_ptr[0].key,"Author");
   text_ptr[1].text = textauthor_;
   text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;
   text_ptr[2].key = new(char[12]); strcpy(text_ptr[0].key,"Description");
   text_ptr[2].text = textdescription_;
   text_ptr[2].compression = PNG_TEXT_COMPRESSION_NONE;
   text_ptr[3].key = new(char[14]); strcpy(text_ptr[0].key,"Creation Time");
   text_ptr[3].text = png_convert_to_rfc1123(png_ptr, &mod_time);
   text_ptr[3].compression = PNG_TEXT_COMPRESSION_NONE;
   text_ptr[4].key = new(char[9]); strcpy(text_ptr[0].key,"Software");
   text_ptr[4].text = textsoftware_;
   text_ptr[4].compression = PNG_TEXT_COMPRESSION_NONE;
   png_set_text(png_ptr, info_ptr, text_ptr, 5);
   
   png_write_info(png_ptr, info_ptr);                     
   png_write_image(png_ptr, graph_);                       
   png_write_end(png_ptr, info_ptr);                      
   png_destroy_write_struct(&png_ptr, &info_ptr);          
   fclose(fp); 
   delete[] text_ptr[0].key;
   delete[] text_ptr[1].key;
   delete[] text_ptr[2].key;
   delete[] text_ptr[3].key;
   delete[] text_ptr[4].key;
}

//////////////////////////////////////////////////////
void pngwriter::line(int xfrom, int yfrom, int xto, int yto, int red, int green,int  blue)
{
   int temp;
   if(xfrom > xto) //Line is always drawn from left to right
     {
	temp = xfrom;
	xfrom= xto;
	xto = temp;
	
	temp = yfrom;      //And if we flipped it in x, must flip it in y too.
	yfrom= yto;
	yto = temp;
	
     }
   
   if(xfrom == xto)
     {
	verticalline(xfrom,yfrom,yto,red,green,blue);
	return;
     }
   
   if(yfrom == yto)
     {
	horizontalline(yfrom,xfrom,xto,red,green,blue);
	return;
     }
   
   double m = ((double)(yto - yfrom))/((double)(xto - xfrom));
   
   if(fabs(m)<1.0)
     {
	yxline(xfrom,yfrom,xto,yto,red,green,blue,m);
	return;
     }
   
   if(fabs(m)>1.0)
     {
	xyline(xfrom,yfrom,xto,yto,red,green,blue,m);
	return;
     }   
   
   if((m==-1.0)&&(xfrom<xto))
     {
	for(temp = 0;temp < fabs(double(xfrom-xto))+1;temp++)
	  {
	     plot(xfrom+temp,yfrom-temp,red,green,blue);
	  }
	return;
     }
   
   if((m==-1.0)&&(xfrom>xto))
     {
	for(temp = 0;temp < fabs(double(xfrom-xto))+1;temp++)
	  {
	     plot(xto+temp,yto-temp,red,green,blue);
	  }
	return;
     }
   
	
   
   if(m==1.0)
     {
	for(temp = 0;temp < fabs(double(xfrom-xto))+1;temp++)
	  {
	     plot(xfrom+temp,yfrom+temp,red,green,blue);
	  }
	return;
     }
   
}

void pngwriter::line(int xfrom, int yfrom, int xto, int yto, double red, double green,double  blue)
{ 
   pngwriter::line( xfrom, 
	 yfrom,
	 xto,
	 yto,
	int (red*65535),
	int (green*65535),
	int (blue*65535)
	);
}


////////////////////////////////////////////////////////////////////////////////////////////
void pngwriter::horizontalline(int yfrom, int xfrom, int xto, int red, int green, int blue)
{
   for(int i=xfrom; i<xto+1; i++)
     {
	plot(i,yfrom,red,green,blue);
     }
}


////////////////////////////////////////////////////////////////////////////////////////////
void pngwriter::verticalline(int xfrom, int yfrom, int yto, int red, int green, int blue)
{
   if(yfrom<yto)
     {
	for(int ii=yfrom; ii<yto+1; ii++)
	  {
	     plot(xfrom,ii,red,green,blue);
	  }
     }
   if(yfrom>yto)
     {
	for(int ii=yto; ii<yfrom+1; ii++)
	  {
	     plot(xfrom,ii,red,green,blue);
	  }
     }
   else
     {
	plot(xfrom,yfrom,red,green,blue);
     }
   
}


////////////////////////////////////////////////////////////////////////////////////////////
void pngwriter::yxline(int xfrom1, int yfrom1, int xto1, int yto1, int red1, int green1, int blue1, double mmm)
{
   int y;

   for(int ip=xfrom1; ip < xto1+1 ; ip++)
     {
	y=((int) (yfrom1 + mmm*((double)(ip-xfrom1))));
	plot(ip,(int)y,red1,green1,blue1);
     }
   
}


////////////////////////////////////////////////////////////////////////////////////////////
void pngwriter::xyline(int xfrom1, int yfrom1, int xto1, int yto1, int red1, int green1, int blue1, double mmm)
{
   
   int x;

   if(yto1>yfrom1)
     {
	for(int ki=yfrom1;  ki < yto1 + 1;ki++)
	  {
	     x=((int)(xfrom1 + (1/mmm)*((double)(ki-yfrom1))));
	     plot(x,ki,red1,green1,blue1);
	  }
     }
   
   if(yto1<yfrom1)
     {
	for(int ki=yto1;  ki < yfrom1 + 1;ki++)
	  {
	     x=((int)(xto1 + (1/mmm)*((double)(ki-yto1))));
	     plot(x,ki,red1,green1,blue1);
	  }
     }	
}

///////////////////////////////////////////////////////////////////////////////////////////
void pngwriter::square(int xfrom, int yfrom, int xto, int yto, int red, int green, int blue)
{
   line(xfrom, yfrom, xfrom, yto, red, green, blue);
   line(xto, yfrom, xto, yto, red, green, blue);
   line(xfrom, yfrom, xto, yfrom, red, green, blue);
   line(xfrom, yto, xto, yto, red, green, blue);      
}

void pngwriter::square(int xfrom, int yfrom, int xto, int yto, double red, double green, double blue)
{
   square( xfrom,  yfrom,  xto,  yto, int(red*65535), int(green*65535), int(blue*65535));
}

//////////////////////////////////////////////////////////////////////////////////////////////////
void pngwriter::filledsquare(int xfrom, int yfrom, int xto, int yto, int red, int green, int blue)
{
   for(int caca = xfrom; caca <xto+1; caca++)
     {
	line(caca, yfrom, caca, yto, red, green, blue);
     }
}

void pngwriter::filledsquare(int xfrom, int yfrom, int xto, int yto, double red, double green, double blue)
{
   filledsquare( xfrom,  yfrom,  xto,  yto, int(red*65535), int(green*65535), int(blue*65535));
}


//////////////////////////////////////////////////////////////////////////////////////////////////   
void pngwriter::circle(int xcentre, int ycentre, int radius, int red, int green, int blue)
{
   if(havesinecosinetables_ == 0)
     {
	if(width_>height_)
	  {
	     sinetable_ = new double[10*width_];
	     cosinetable_ = new double[10*width_];
	     
	     for(int  uuu = 0; uuu < 10*width_; uuu++)
	       {
		  sinetable_[uuu] = sin( (M_PI+M_PI)*((double)uuu)/(10.0*((double)width_)));
		  cosinetable_[uuu] = cos( (M_PI+M_PI)*((double)uuu)/(10.0*((double)width_)));
	       }
	  }
	else
	  {
	     sinetable_ = new double[10*height_];
	     cosinetable_ = new double[10*height_];
	     
	     for(int  uuu = 0; uuu < 10*height_; uuu++)
	       {
		  sinetable_[uuu] = sin( (M_PI+M_PI)*((double)uuu)/(10.0*((double)height_)));
		  cosinetable_[uuu] = cos( (M_PI+M_PI)*((double)uuu)/(10.0*((double)height_)));
	       }
	  }
     }//if havesinecosinetables_ == 0
   
   if(width_>height_)
     {
	for(int jjj = 0; jjj< 10*width_; jjj++)
	  {
	     plot(xcentre + int(sinetable_[jjj]*((double)radius)),ycentre + int(cosinetable_[jjj]*((double)radius)),red,green,blue);
	  }
     }
   else
       {
	for(int jjj = 0; jjj< 10*height_; jjj++)
	  {
	     plot( ((int)(xcentre + sinetable_[jjj]*((double)radius))),((int)(ycentre + cosinetable_[jjj]*((double)radius))),red,green,blue);
	  }
     }
   
}

void pngwriter::circle(int xcentre, int ycentre, int radius, double red, double green, double blue)
{
   circle(xcentre,ycentre,radius, int(red*65535), int(green*65535), int(blue*65535));   
}

////////////////////////////////////////////////////////////
void pngwriter::filledcircle(int xcentre, int ycentre, int radius, int red, int green, int blue)
{
   for(int jjj = ycentre-radius; jjj< ycentre+radius+1; jjj++)
     {
	line(xcentre - int(sqrt((double)(radius*radius) - (-ycentre + jjj)*(-ycentre + jjj ))), jjj, 
	     xcentre + int(sqrt((double)(radius*radius) - (-ycentre + jjj)*(-ycentre + jjj ))),jjj,red,green,blue);
     }
}

void pngwriter::filledcircle(int xcentre, int ycentre, int radius, double red, double green, double blue)
{
  filledcircle( xcentre, ycentre,  radius, int(red*65535), int(green*65535), int(blue*65535));
}

////////////////Reading routines/////////////////////
/////////////////////////////////////////////////

void pngwriter::readfromfile(char * name)
{
   FILE            *fp;
   png_structp     png_ptr;
   png_infop       info_ptr;
   unsigned char   **image;
   unsigned long   width, height;
   int bit_depth, color_type, interlace_type;
   //   png_uint_32     i;
   
   check_if_png(name, &fp);
   read_png_info(fp, &png_ptr, &info_ptr);
   read_png_image(fp, png_ptr, info_ptr, &image, &width, &height);
   //stuf should now be in image[][].
   
   //First we must get rid of the image already there, and free the memory.
   int jjj;
   for (jjj = 0; jjj < height_; jjj++) free(graph_[jjj]);            
   free(graph_);
      
   //Must reassign the new size of the read image
   width_ = width;
   height_ = height;
   
   //Graph now is the image.
   graph_ = image;
   
   rowbytes_ = png_get_rowbytes(png_ptr, info_ptr);
   
   png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type, &interlace_type, NULL, NULL);
   bit_depth_ = bit_depth;
   colortype_ = color_type;

   if(color_type == PNG_COLOR_TYPE_PALETTE /*&& bit_depth<8*/)
     {
	png_set_expand(png_ptr);
     }
   
   if(color_type == PNG_COLOR_TYPE_GRAY && bit_depth<8)
     {
	png_set_expand(png_ptr);
     }
   
   if(color_type & PNG_COLOR_MASK_ALPHA)
     {
	png_set_strip_alpha(png_ptr);
     }
   
   if(color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_RGB_ALPHA)
     {
	png_set_gray_to_rgb(png_ptr);
     }
   
   if((bit_depth_ !=16)&&(bit_depth_ !=8))
     {
	cerr << "PNGwriter::readfromfile() - WARNING **: Input file is of unsupported type (bad bit_depth). Output will be unpredictable.\n";
	cout << bit_depth_ << endl;
     }
  
    if(colortype_ !=2)
     {
	cerr << "PNGwriter::readfromfile() - WARNING **: Input file is of unsupported type (bad color_type). Output will be unpredictable.\n";
     }
  
     
   screengamma_ = 2.2;  
   double          file_gamma,screen_gamma;
   screen_gamma = screengamma_;
   if (png_get_gAMA(png_ptr, info_ptr, &file_gamma))
     {
//	cout << "filegamma from info is " << file_gamma << endl;
	png_set_gamma(png_ptr,screen_gamma,file_gamma);
     }
   else
     {
	png_set_gamma(png_ptr, screen_gamma,0.45);
     }
   
   
   
   filegamma_ = file_gamma; 
}

///////////////////////////////////////////////////////

void pngwriter::readfromfile(const char * name)
{
   readfromfile((char *)(name));
}




/////////////////////////////////////////////////////////
void pngwriter::check_if_png(char *file_name, FILE **fp)
{
   char    sig[PNG_BYTES_TO_CHECK];
   
   if ((*fp = fopen(file_name, "rb")) == NULL)
     {
	//   exit(EXIT_FAILURE);
   	cerr << "PNGwriter::check_if_png - ERROR **: Could not open file to read." << endl;
     }
   
   if (fread(sig, 1, PNG_BYTES_TO_CHECK, *fp) != PNG_BYTES_TO_CHECK) 
     {
	fclose(*fp);
	//exit(EXIT_FAILURE);
	cerr << "PNGwriter::check_if_png - ERROR **: File does not appear to be a valid PNG file." << endl;
     }
}


///////////////////////////////////////////////////////
void pngwriter::read_png_info(FILE *fp, png_structp *png_ptr, png_infop *info_ptr)
{
   *png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
   if (*png_ptr == NULL) 
     {
	fclose(fp);
	cerr << "PNGwriter::read_png_info - ERROR **: Could not create read_struct." << endl;
	//exit(EXIT_FAILURE);
     }
   *info_ptr = png_create_info_struct(*png_ptr);
   if (*info_ptr == NULL) 
     {
	png_destroy_read_struct(png_ptr, (png_infopp)NULL, (png_infopp)NULL);
	fclose(fp);
	cerr << "PNGwriter::read_png_info - ERROR **: Could not create info_struct." << endl;
	//exit(EXIT_FAILURE);
     }
   if (setjmp((*png_ptr)->jmpbuf)) 
     {
	png_destroy_read_struct(png_ptr, info_ptr, (png_infopp)NULL);
	fclose(fp);
	cerr << "PNGwriter::read_png_info - ERROR **: Unga-bunga not the right shade of pink (setjmp(*png_ptr)->jmpbf) failed)." << endl;
	//exit(EXIT_FAILURE);
     }
   png_init_io(*png_ptr, fp);
   png_set_sig_bytes(*png_ptr, PNG_BYTES_TO_CHECK);
   png_read_info(*png_ptr, *info_ptr);
}


////////////////////////////////////////////////////////////
void pngwriter::read_png_image(FILE *fp, png_structp png_ptr, png_infop info_ptr,
		    png_bytepp *image, png_uint_32 *width, png_uint_32 *height)
{
   unsigned int i,j;
   
   *width = png_get_image_width(png_ptr, info_ptr);
   *height = png_get_image_height(png_ptr, info_ptr);
   if ((*image = (png_bytepp)malloc(*height * sizeof(png_bytep))) == NULL) 
     {
	fclose(fp);
	cerr << "PNGwriter::read_png_image - ERROR **: Could not allocate memory for reading image." << endl;
	//exit(EXIT_FAILURE);
     }
   for (i = 0; i < *height; i++)
     {
	(*image)[i] = (png_bytep)malloc(png_get_rowbytes(png_ptr, info_ptr));
	if ((*image)[i] == NULL)
	  {
	     for (j = 0; j < i; j++) free((*image)[j]);
	     free(*image);
	     fclose(fp);
	     cerr << "PNGwriter::read_png_image - ERROR **: Could not allocate memory for reading image." << endl;
	     //exit(EXIT_FAILURE);
	  }
     }
   png_read_image(png_ptr, *image);
}


///////////////////////////////////
int pngwriter::getheight(void)
{
   return height_;
}


int pngwriter::getwidth(void)
{
   return width_;
}

int pngwriter::getbitdepth(void)
{
   return bit_depth_;
}


int pngwriter::getcolortype(void)
{
   return colortype_;
}


double pngwriter::getgamma(void)
{
   return filegamma_;
}


void pngwriter::setgamma(double gamma)
{
    filegamma_ = gamma;
}

// The algorithms HSVtoRGB and RGBtoHSV were found at http://www.cs.rit.edu/~ncs/
//  which is a page that belongs to Nan C. Schaller, though
//  these algorithms appear to be the work of Eugene Vishnevsky. 
//////////////////////////////////////////////
void pngwriter::HSVtoRGB( double *r, double *g, double *b, double h, double s, double v )
{	     
   // r,g,b values are from 0 to 1 
   // h = [0,1], s = [0,1], v = [0,1] 
   // if s == 0, then h = -1 (undefined)
  
   h = h*360.0;
   
   int i; 
   double f, p, q, t;
    if( s == 0 ) 
     {
	// achromatic (grey) 
	*r = *g = *b = v; 
	return; 
     }
   
   h /= 60;                        // sector 0 to 5 
   i = int(floor( h )); 
   f = h - i;                      // factorial part of h 
   p = v * ( 1 - s ); 
   q = v * ( 1 - s * f ); 
   t = v * ( 1 - s * ( 1 - f ) );
   
   switch( i ) { 
    case 0: 
      *r = v; 
      *g = t; 
      *b = p; 
      break; 
    case 1: 
      *r = q; 
      *g = v; 
      *b = p; 
      break; 
    case 2: 
	*r = p; 
      *g = v; 
      *b = t; 
      break; 
    case 3: 
      *r = p; 
      *g = q; 
      *b = v; 
      break; 
    case 4: 
      *r = t; 
      *g = p; 
      *b = v; 
      break; 
    default:                // case 5: 
      *r = v; 
      *g = p; 
      *b = q; 
      break; 
   }   
}

void pngwriter::RGBtoHSV( float r, float g, float b, float *h, float *s, float *v ) 
{
   
   float min=0.0, max=0.0, delta;
   

   if( (r>=g)&&(r>=b) )
       {
	  max = r;
       }
   if( (g>=r)&&(g>=b) )
       {
	  max = g;
       }
   if( (b>=g)&&(b>=r) )
       {
	  max = b;
       }
   
   if( (r<=g)&&(r<=b) )
       {
	  min = r;
       }
   if( (g<=r)&&(g<=b) )
       {
	  min = g;
       }
   if( (b<=g)&&(b<=r) )
       {
	  min = b;
       }   
   
   *v = max;                               // v
   
   delta = max - min;
   
   if( max != 0 ) 
     *s = delta / max;               // s 
   else 
     {
	
	r = g = b = 0;                // s = 0, v is undefined 
	  *s = 0; 
	*h = -1; 
	return; 
     }
   
   if( r == max ) 
     *h = ( g - b ) / delta;         // between yellow & magenta 
   else if( g == max ) 
     *h = 2 + ( b - r ) / delta;     // between cyan & yellow 
   else 
     *h = 4 + ( r - g ) / delta;     // between magenta & cyan
   
   *h *= 60;                               // degrees 
   if( *h < 0 ) 
     *h += 360;
   
}


//
//////////////////////////////////////////////////////////////////////////////////
void pngwriter::plotHSV(int x, int y, double hue, double saturation, double value)
{
   double red,green,blue;
   double *redp;
   double *greenp;
   double *bluep;
    
   redp = &red;
   greenp = &green;
   bluep = &blue;
 
   HSVtoRGB(redp,greenp,bluep,hue,saturation,value);   
   plot(x,y,red,green,blue);
}

void pngwriter::plotHSV(int x, int y, int hue, int saturation, int value)
{
   plotHSV(x, y, double(hue)/65535.0, double(saturation)/65535.0,  double(value)/65535.0);
}



//
//////////////////////////////////////////////////////////////////////////////////
double pngwriter::dreadHSV(int x, int y, int colour)
{  
   float * huep;
   float * saturationp;
   float * valuep;
   float red,green,blue;
   float hue, saturation, value;
   
   red = float(dread(x,y,1));
   green = float(dread(x,y,2));
   blue = float(dread(x,y,3));
   
   huep = &hue;
   saturationp = &saturation;
   valuep = &value;
 
   RGBtoHSV( red,  green,  blue, huep,  saturationp, valuep );
   
   if(colour == 1)
     {
	return double(hue)/360.0;
     }
   
   else if(colour == 2)
     {
	return saturation;
     }
   
   else if(colour == 3)
     {
	return value;
     }
   else
     {
	cerr << "PNGwriter::dreadHSV - ERROR **: Called with wrong colour argument: should be 1, 2 or 3; was: " << colour << "." << endl;   
     }

   return 0; // should never get here
}


//
//////////////////////////////////////////////////////////////////////////////////
int pngwriter::readHSV(int x, int y, int colour)
{
   float * huep;
   float * saturationp;
   float * valuep;
   float red,green,blue;
   float hue, saturation, value;
   
   red = float(dread(x,y,1));
   green = float(dread(x,y,2));
   blue = float(dread(x,y,3));
   
   huep = &hue;
   saturationp = &saturation;
   valuep = &value;
 
   RGBtoHSV( red,  green,  blue, huep,  saturationp, valuep );
   
   if(colour == 1)
     {
	return int(65535*(double(hue)/360.0));
     }
   
   else if(colour == 2)
     {
	return int(65535*saturation);
     }
   
   else if(colour == 3)
     {
	return int(65535*value);
     }
   else
     {
	cerr << "PNGwriter::readHSV - ERROR **: Called with wrong colour argument: should be 1, 2 or 3; was: " << colour << "." << endl;   
     }

   return 0; // should never get here
}

void pngwriter::setcompressionlevel(int level)
{
   if( (level < -1)||(level > 9) )
     {
	cerr << "PNGwriter::setcompressionlevel - ERROR **: Called with wrong compression level: should be -1 to 9, was: " << level << "." << endl;   
     }
   compressionlevel_ = level;  
}


// An implementation of a Bezier curve.
void pngwriter::bezier(  int startPtX, int startPtY,                                                                            
              int startControlX, int startControlY,                                                                                     
              int endPtX, int endPtY,                                                                                                   
              int endControlX, int endControlY,                                                                                         
              double red, double green, double blue)                                                                                    
{
                                                                                                                                        
      double cx = 3.0*(startControlX - startPtX);                                                                                        
      double bx = 3.0*(endControlX - startControlX) - cx;                                                                                
      double ax = double(endPtX - startPtX - cx - bx);                                                                                   
                                                                                                                                         
      double cy = 3.0*(startControlY - startPtY);                                                                                        
      double by = 3.0*(endControlY - startControlY) - cy;                                                                                
      double ay = double(endPtY - startPtY - cy - by);                                                                                   
                                                                                                                                         
      double x,y,newx,newy;                                                                                                              
      x = startPtX;                                                                                                                      
      y = startPtY;                                                                                                                      
                                                                                                                                         
      for(double t = 0.0; t<=1.005; t += 0.005)                                                                                          
     {
	newx = startPtX + t*(double(cx) + t*(double(bx) + t*(double(ax))));                                                           
	newy = startPtY + t*(double(cy) + t*(double(by) + t*(double(ay))));                                                           
	this->line(int(x),int(y),int(newx),int(newy),red,green,blue);                                                                
	x = newx;                                                                                                                     
	y = newy;                                                                                                                     
     }                                                                                                                                   
}
      
//int version of bezier
void pngwriter::bezier(  int startPtX, int startPtY,                                                                            
              int startControlX, int startControlY,                                                                                     
              int endPtX, int endPtY,                                                                                                   
              int endControlX, int endControlY,                                                                                         
              int  red, int  green, int blue)                                                                                    
{
 bezier(   startPtX,  startPtY,                                                                            
              startControlX, startControlY,                                                                                     
              endPtX, endPtY,                                                                                                   
              endControlX,  endControlY, 
              double(red)/65535.0,  double(green)/65535.0,  double(blue)/65535.0); 
}

/*
int pngwriter::getcompressionlevel(void)
{
   return png_get_compression_level(png_ptr);
}
*/
