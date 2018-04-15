/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
/*
Simple progressbar modified from a code by Remik Ziemlinski
See below for more infos.
Nguyen Damien 2013

Copyright (C) 2011,2012 Remik Ziemlinski. See MIT-LICENSE.

CHANGELOG

v0.0.0 20110502 rsz Created.
V1.0.0 20110522 rsz Extended to show eta with growing bar.
v2.0.0 20110525 rsz Added time elapsed.
v2.0.1 20111006 rsz Added default constructor value.
*/

#ifndef PROGRESS_BAR_HPP_INCLUDED
#define PROGRESS_BAR_HPP_INCLUDED

#include <iostream>
#include <string>
#include <stdio.h>

#pragma GCC system_header


#ifdef SIMPLE_PROGRESS_BAR
// Very simple, one line (without refreshing) progress bar
// 0...10...20...30...40...50...60...70...80...90...100
class Progress_Bar {
public:
     Progress_Bar(unsigned int _n=0, std::string prefix_a = "") 
	  : n(_n), pct(0), cur(0), prefix(prefix_a)
	  {}
     void reset(unsigned int n_a = 0,
		std::string prefix_a = "") { 
	       n = n_a; pct = 0; cur = 0; prefix = prefix_a;
	  }

     void start() { 
	  std::cout << prefix << '0'; std::cout.flush(); 
     }

     void operator++() {
	  if (cur >= n) {
	       return;
	  }
	  ++cur;
	  setPct( (float)cur/n );
     };
	
     // Set 0.0-1.0, where 1.0 equals 100%.
     void setPct(float Pct) {
	  short delta = (short)(Pct*1000 - pct);
	  if (delta < 25) return;
		
	  do {
	       pct += 25;
	       if ( (pct % 100) == 0 ) 
		    std::cout << pct/10;
	       else
		    std::cout << '.';
	  } while((delta -= 25) >= 25);

	  if (Pct >= 1.0) {
	       std::cout << std::endl;
	  }

	  std::cout.flush();
     };

     unsigned int n;
     unsigned int cur;
     unsigned short pct; // Stored as 0-1000, so 2.5% is encoded as 25.
     std::string prefix;
};

#else 
// One-line refreshing progress bar inspired by wget
// 90% [##################################################     ]
class Progress_Bar {
public:
     Progress_Bar(unsigned int n_a = 0,
		  std::string prefix_a = "")
	  : n(n_a), pct(0), cur(0), width(80),
	    prefix(prefix_a)
     {}
	  
     void reset(unsigned int n_a = 0,
		std::string prefix_a = "") { 
	  n = n_a; pct = 0; cur = 0; prefix = prefix_a;}
     void start() { setPct(0); }
	
     void operator++() {
	  if (cur >= n) return;
	  ++cur;
		
	  setPct( ((float)cur)/n );
     }

     // Set 0.0-1.0, where 1.0 equals 100%.
     void setPct(float Pct) {
	  char pctstr[5];
	  std::cout << '\r' << prefix;
	  sprintf(pctstr, "%3d%%", (int)(100*Pct));
	  // Compute how many tics we can display.
	  int nticsMax = (width-27);
	  int ntics = (int)(nticsMax*Pct);
	  std::string out(pctstr);
	  out.append(" [");
	  out.append(ntics,'#');
	  out.append(nticsMax-ntics,' ');
	  out.append("] ");
	  std::string tstr;
	  if (Pct >= 1.0) {
	       out.append("\n");
	       std::cout << out;
	       return;
	  } 

	  // Pad end with spaces to overwrite previous string that may have been longer.
	  if (out.size() < width)
	       out.append(width-out.size(),' ');
			
	  std::cout << out;
	  std::cout.flush();
     }

     unsigned int n;
     unsigned short pct; // Stored as 0-1000, so 2.5% is encoded as 25.
     unsigned int cur;
     unsigned char width; // How many chars the entire line can be.
     std::string prefix;
};
#endif // SIMPLE_PROGRESS_BAR

#endif //PROGRESS_BAR_HPP_INCLUDED
