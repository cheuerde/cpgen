#ifndef _PRINTER_H
#define _PRINTER_H
#include <stdio.h>
#include <Rcpp.h>


// taken from: http://www.codeproject.com/Tips/537904/Console-simple-progress

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


	public:
		inline void DoProgress() {

			if(!initialized) { initialize();}

			progress++;

			percent = ( progress * 100 ) / total;

			if(percent>100) { percent=100; }

			if(percent > previous_percent) {

				fflush(stdout);
				previous_pos = pos;
				pos = ( progress * width ) / total;  
				int step = pos - previous_pos;

				for ( int i = 0; i < step; i++ )  Rcpp::Rcout << "=";

				for(int i=0;i<(width-pos+2);i++) { Rcpp::Rcout << "\033[C"; }

				Rcpp::Rcout << percent << "%" << "\033[" << pos+2 << "G";

				previous_percent = percent;

				if(pos == width) { Rcpp::Rcout << std::endl << std::endl; }

			}


		};

		void initialize() {



			Rcpp::Rcout << std::endl << "[";

			//fill progress bar 
			for ( int i = 0; i < width; i++ )  Rcpp::Rcout << " ";

			Rcpp::Rcout << "]" << " " << "0%"<< "\033[2G" <<  std::flush;

			pos = previous_pos = percent = previous_percent = progress = 0;

			initialized = true;

		};


		printer(int t) : total(t) {} 



};



#endif
