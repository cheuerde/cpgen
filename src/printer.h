/*
// printer.h
// Claas Heuer, June 2014
//
// Copyright (C)  2014 Claas Heuer
//
// This file is part of cpgen.
//
// cpgen is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// cpgen is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// in file.path(R.home("share"), "licenses").  If not, see
// <http://www.gnu.org/licenses/>.
*/

//
// this is a class for printing progress to the screen.
// an option was RcppProgress but it needs some adjustment in terms
// of flushing (unbuffered output).
// Here we simply construct a message using a std::string, jump to the beginning
// of always the same line and print the message. A more elegant
// solution with ANSI-escapes is in 'printer_unix.h'. Unfortunately,
// ANSI-escapes are not supported by Windows, that is why we print the whole
// line over and over again, because the only thing that works in unix and windows
// is the line-escape '\r'.
//
// Main idea was taken from: http://www.codeproject.com/Tips/537904/Console-simple-progress
//

#ifndef _PRINTER_H
#define _PRINTER_H
#include <stdio.h>
#include <iomanip>
#include <Rcpp.h>
#include <iostream> 
#include <string.h>



class printer {

	private:
		int percent;
		int previous_percent;
		static const int width = 40;
		int pos;
		int previous_pos;
		int total;
		bool initialized;
		int progress;
		std::string message;


	public:
		inline void DoProgress() {

			if(!initialized) { initialize();}

			progress++;

			percent = ( progress * 100 ) / total;

			if(percent>100) { percent=100; }

			if(percent > previous_percent) {


				previous_pos = pos;
				pos = ( progress * width ) / total;  
				//      int step = pos - previous_pos;

				std::ostringstream oss;
				oss << percent;

				message=" [";
				for(int i=0;i<pos;i++) { message.append("="); }
				for(int i=0;i<width-pos;i++)  {message.append(" "); }
				message.append("] ");
				message.append(oss.str());
				message.append("%");

				Rcpp::Rcout << message;
				Rcpp::Rcout << "\r" << std::flush; 
				//      fflush(stdout);

				previous_percent = percent;

				if(pos == width) { 

					Rcpp::Rcout << std::endl; 
					Rcpp::Rcout << "\r  " << std::flush;
					Rcpp::Rcout << std::endl; 

				}

			}


		};

		void initialize() {



			Rcpp::Rcout << std::endl << " [";

			//fill progress bar 
			for ( int i = 0; i < width; i++ )  Rcpp::Rcout << " ";

			Rcpp::Rcout << "]" << " " << "0%"<< "\r" <<  std::flush;

			initialized = true;

		};


		printer(int t) : total(t) {


			initialized=false;
			pos = previous_pos = percent = previous_percent = progress = 0;

		} 



};



#endif
