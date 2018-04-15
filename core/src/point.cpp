#include "point.hpp"

cellsim::Point::Point()
     : Point(0., 0., 0.)
{}

cellsim::Point::Point(double c_in, double I_in, double s_in)
     : c_(c_in), I_(I_in), s_(s_in)
{
     assert(c_ >= 0.);
     assert(I_ >= 0.);
     assert(s_ >= 0.);
}

std::ostream& operator<<(std::ostream& out, cellsim::Point p)
{
     out << p.c() << ", " << p.I() << ", " << p.s();
     return out;
}
