#ifndef AB_EXCEPTION_H
#define AB_EXCEPTION_H
/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

//DESCRIPTION:
//This is a simple but useful class for construction of standard exceptions
//The constructor is implemented to take variadic arguments 
//in the style of printf for simple message construction

#include <exception>
#include <string>
#include <cstdarg>

class MyException : public std::exception{
	public:
		MyException(const std::string &m):msg(m), c(0){}
		//my first variadic constructor :)
		MyException(const char *str, ...):c(0){
			char *buffer = new char[1000];
			va_list ap;
			va_start(ap, str);
			vsprintf(buffer, str, ap);
			va_end(ap);
			msg = buffer;
			delete [] buffer;
		}
		MyException(const int type, const char *str, ...):c(type){
			char *buffer = new char[1000];
			va_list ap;
			va_start(ap, str);
			vsprintf(buffer, str, ap);
			va_end(ap);
			msg = buffer;
			delete [] buffer;
		}
		virtual ~MyException() throw() {};
		virtual const char* what() const throw(){
			return msg.c_str();
		}
		virtual const int code()const{ return c; }
	private:
		std::string msg;	
		int c;
};
#endif
