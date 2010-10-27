//**********  pngwriter.h   **********************************************
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
//  Documentation:             This header file is commented, but for a
//                             quick reference document, and support,
//                             take a look at the website.
//
//*************************************************************************

#ifndef PNGWRITER_H
#define PNGWRITER_H 1

#define PNGWRITER_VERSION 0.36

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <png.h>
#include <math.h>
#include <setjmp.h>
#include <string>
using namespace std;
using std::ios;
using std::ifstream;
using std::ofstream;

#define PNG_BYTES_TO_CHECK (4)
#define PNGWRITER_DEFAULT_COMPRESSION (6)

class pngwriter 
{
 private:
   
   char filename_[255];   
   int height_;
   int width_;
   int  backgroundcolour_;
   int bit_depth_;
   int rowbytes_;
   int colortype_;
   int compressionlevel_;
   unsigned char * * graph_;
   double filegamma_;
   double screengamma_;
   bool havesinecosinetables_;
   double * sinetable_;
   double * cosinetable_;
   char textauthor_[255];   
   char textdescription_[255];   
   char texttitle_[255];   
   char textsoftware_[255];   
   void check_if_png(char *file_name, FILE **fp);
   void read_png_info(FILE *fp, png_structp *png_ptr, png_infop *info_ptr);
   void read_png_image(FILE *fp, png_structp png_ptr, png_infop info_ptr,
 		       png_bytepp *image, png_uint_32 *width, png_uint_32 *height);
   void horizontalline(int yfrom, int xfrom, int xto, int red, int green, int blue);
   void verticalline(int xfrom, int yfrom, int yto, int red, int green, int blue);
   void yxline(int xfrom, int yfrom, int xto, int yto, int red, int green, int blue, double mmm);
   void xyline(int xfrom, int yfrom, int xto, int yto, int red, int green, int blue, double mmm);
   
   /* The algorithms HSVtoRGB and RGBtoHSV were found at http://www.cs.rit.edu/~ncs/
    * which is a page that belongs to Nan C. Schaller, though
    * these algorithms appear to be the work of Eugene Vishnevsky. 
    * */
   void HSVtoRGB( double *r, double *g, double *b, double h, double s, double v ); 
   void RGBtoHSV( float r, float g, float b, float *h, float *s, float *v );
   
 public:

   /* General Notes
    * It is important to remember that all functions that accept an argument of type "const char *" will also
    * accept "char *", this is done so you can have a changing filename (to make many PNG images in series 
    * with a different name, for example), and to allow you to use string type objects which can be easily 
    * turned into const char * (if theString is an object of type string, then it can be used as a const char *
    * by saying theString.c_str()).
    * It is also important to remember that whenever a function has a colour coeffiecient as its argument, 
    * that argument can be either an int from 0 to 65535 or a double from 0.0 to 1.0. 
    * It is important to make sure that you are calling the function with the type that you want.
    * Remember that 1 is an int, while 1.0 is a double, and will thus determine what version of the function 
    * will be used. Similarly, do not make the mistake of calling for example plot(x, y, 0.0, 0.0, 65535),
    * because
    * there is no plot(int, int, double, double, int).
    * */
   
   /* Constructor
    * The constructor requires the width and the height of the image, the background colour for the
    * image and the filename of the file (a pointer or simple "myfile.png"). The background colour
    * can only be initialized to a shade of grey (once the object has been created you can do whatever 
    * you want, though), because generally one wants either a white (65535 or 1.0) or a black (0 or 0.0)
    * background to start with.
    * Tip: The filename can be given as easily as:
    * pngwriter mypng(300, 300, 0.0, "myfile.png");    
    * Tip: If you are going to create a PNGwriter instance for reading in a file that already exists, 
    * then width and height can be 1 pixel, and the size will be automatically adjusted once you use
    * readfromfile().
    * */
   pngwriter(int width, int height, int backgroundcolour, char * filename);   
   pngwriter(int width, int height, double backgroundcolour, char * filename);    
   pngwriter(int width, int height, int backgroundcolour, const char * filename);   
   pngwriter(int width, int height, double backgroundcolour, const char * filename);    

   /* Destructor
    * */
   ~pngwriter();  
   
   /*  Plot
    * With this function a pixel at coordinates (x, y) can be set to the desired colour. 
    * The pixels are numbered starting from (1, 1) and go to (width, height). 
    * As with most functions in PNGwriter, it has been overloaded to accept either int arguments 
    * for the colour coefficients, or those of type double. If they are of type int, 
    * they go from 0 to 65535. If they are of type double, they go from 0.0 to 1.0.
    * Tip: To plot using red, then specify plot(x, y, 0.0, 0.0, 1.0). To make pink, 
    * just add a constant value to all three coefficients, like this:
    * plot(x, y, 0.4, 0.4, 1.0). 
    * Tip: If nothing is being plotted to your PNG file, make sure that you remember
    * to close() the instance before your program is finished, and that the x and y position
    * is actually within the bounds of your image. If either is not, then PNGwriter will 
    * not complain-- it is up to you to check for this!
    * */
   void  plot(int x, int y, int red, int green, int blue); 
   void  plot(int x, int y, double red, double green, double blue); 
                                                          
   /* Plot HSV
    * With this function a pixel at coordinates (x, y) can be set to the desired colour, 
    * but with the colour coefficients given in the Hue, Saturation, Value colourspace. 
    * This has the advantage that one can determine the colour that will be plotted with 
    * only one parameter, the Hue. The colour coefficients must go from 0 to 65535 and
    * be of type int, or be of type double and go from 0.0 to 1.0.
    * */
   void plotHSV(int x, int y, double hue, double saturation, double value);
   void plotHSV(int x, int y, int hue, int saturation, int value); 

   /* Read
    * With this function we find out what colour the pixel (x, y) is. If "colour" is 1,
    * it will return the red coefficient, if it is set to 2, the green one, and if 
    * it set to 3, the blue colour coefficient will be returned,
    * and this returned value will be of type int and be between 0 and 65535.
    * */
   int read(int x, int y, int colour);  

   /* Read, Average
    * Same as the above, only that the average of the three colour coefficients is returned.
    */
   int read(int x, int y);              
   
  /* dRead
   * With this function we find out what colour the pixel (x, y) is. If "colour" is 1,
   * it will return the red coefficient, if it is set to 2, the green one, and if 
   * it set to 3, the blue colour coefficient will be returned,
   * and this returned value will be of type double and be between 0.0 and 1.0.
   * */
   double dread(int x, int y, int colour);   
   
   /* dRead, Average
    * Same as the above, only that the average of the three colour coefficients is returned.
    */
   double dread(int x, int y);             
   
   /* Read HSV
    * With this function we find out what colour the pixel (x, y) is, but in the Hue, 
    * Saturation, Value colourspace. If "colour" is 1,
    * it will return the Hue coefficient, if it is set to 2, the Saturation one, and if 
    * it set to 3, the Value colour coefficient will be returned,
    * and this returned value will be of type int and be between 0 and 65535.
    * Tip: This is especially useful for categorizing sections of the image according 
    * to their colour. 
    * */
   int readHSV(int x, int y, int colour);  
    
  /* dRead HSV
   * With this function we find out what colour the pixel (x, y) is, but in the Hue, 
   * Saturation, Value colourspace. If "colour" is 1,
    * it will return the Hue coefficient, if it is set to 2, the Saturation one, and if 
    * it set to 3, the Value colour coefficient will be returned,
    * and this returned value will be of type double and be between 0.0 and 1.0.
    * */
   double dreadHSV(int x, int y, int colour);    

   /* Clear
    * The whole image is set to black.
    * */ 
   void clear(void);    
   
   /* Close
    * Close the instance of the class, and write the image to disk.
    * Tip: If you do not call this function before your program ends, no image
    * will be written to disk.
    * */
   void close(void); 

   /* Rename
    * To rename the file once an instance of pngwriter has been created.
    * Useful for assigning names to files based upon their content.
    * Tip: This is as easy as calling pngwriter_rename("newname.png")
    * */
   void pngwriter_rename(char * newname);               
   void pngwriter_rename(const char * newname);            

   /* Figures
    * These functions draw basic shapes. Available in both int and double version.
    * They are fast and efficient. For example, there are several ways to program
    * the drawing of a line. Some are efficient computationally, and some are not.
    * The one that line uses is the most efficient. Another example is how the functions
    * that draw circles are implemented. Calling sin() and cos() is very expensive
    * computationally, so the first time that these are called, PNGwriter generates
    * an internal table of sines and cosines. For all subsequent calls to circle() 
    * or filledcircle() the values of the table will be used.
    * */
   void line(int xfrom, int yfrom, int xto, int yto, int red, int green,int  blue);
   void line(int xfrom, int yfrom, int xto, int yto, double red, double green,double  blue);
  
   void square(int xfrom, int yfrom, int xto, int yto, int red, int green,int  blue);
   void square(int xfrom, int yfrom, int xto, int yto, double red, double green,double  blue);
   
   void filledsquare(int xfrom, int yfrom, int xto, int yto, int red, int green,int  blue);
   void filledsquare(int xfrom, int yfrom, int xto, int yto, double red, double green,double  blue);

   void circle(int xcentre, int ycentre, int radius, int red, int green, int blue);
   void circle(int xcentre, int ycentre, int radius, double red, double green, double blue);
   
   void filledcircle(int xcentre, int ycentre, int radius, int red, int green, int blue);
   void filledcircle(int xcentre, int ycentre, int radius, double red, double green, double blue);

   /* Read From File
    * Open the existing PNG image, and copy it into this instance of the class. Now you can access
    * that image's pixels with read(x, y, colour), you can change the image, or whatever you want.
    * It can be called like this, also: readfromfile("image.png"). It is important to mention that 
    * not all colour types and bit depths are supported. Try and make sure that your PNG image is
    * of bit depth 8 or 16.
    * */
   void readfromfile(char * name);  
   void readfromfile(const char * name); 

   /* Get Height
    * When you open a PNG with readfromfile() you can find out its height with this function.
    * */
   int getheight(void);
   
   /* Get Width
    * When you open a PNG with readfromfile() you can find out its width with this function.
    * */
   int getwidth(void);

   /* Set Compression Level
    * Set the compression level that will be used for the image. -1 is default, 0 is none, 9 is best compression. 
    *  Remember that this will affect how long it will take to close() the image.
    * */
    void setcompressionlevel(int level);

   /* Get Bit Depth
    * When you open a PNG with readfromfile() you can find out its bit depth with this function.
    * Mostly for troubleshooting uses.
    * */
   int getbitdepth(void);
   
   /* Get Colour Type
    * When you open a PNG with readfromfile() you can find out its colour type (libpng categorizes 
    * different styles of image data with this number).
    * Mostly for troubleshooting uses.
    * */
   int getcolortype(void);
   
   /* Set Gamma Coeff
    * Set the image's gamma (file gamma) coefficient. This is experimental, but use it if your image's colours seem too bright
    * or too dark. The default value of 0.5 should be fine. 
    * */
   void setgamma(double gamma);

   
   /* Get Gamma Coeff
    * Get the image's gamma coefficient. This is experimental.
    * */
   double getgamma(void);

   /* Bezier Curve
    * (After Frenchman Pierre BŽzier from Regie Renault)
    * A collection of formulae for describing curved lines 
    * and surfaces, first used in 1972 to model automobile surfaces. 
    *             (from the The Free On-line Dictionary of Computing)
    * See http://www.moshplant.com/direct-or/bezier/ for one of many
    * available descriptions of bezier curves.
    * There are four points used to define the curve: the two endpoints
    * of the curve are called the anchor points, while the other points,
    * which define the actual curvature, are called handles or control points.
    * Moving the handles lets you modify the shape of the curve. 
    * */ 
     
   void bezier(  int startPtX, int startPtY,                                                                            
              int startControlX, int startControlY,                                                                                     
              int endPtX, int endPtY,                                                                                                   
              int endControlX, int endControlY,                                                                                         
              double red, double green, double blue);

   void bezier(  int startPtX, int startPtY,                                                                            
              int startControlX, int startControlY,                                                                                     
              int endPtX, int endPtY,                                                                                                   
              int endControlX, int endControlY,                                                                                         
              int red, int green, int blue);

   /* Set Text
    * Sets the text information in the PNG header. If it is not called, the default is used.
    */
   void settext(char * title, char * author, char * description, char * software);
   void settext(const char * title, const char * author, const char * description, const char * software);

   
};


#endif

