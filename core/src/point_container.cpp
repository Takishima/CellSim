#include "point_container.hpp"

typedef cellsim::Point_Container Point_Container;

Point_Container::Point_Container(size_type n_point)
     : Point_Container(n_point, value_type(), value_type(), value_type())
{}

Point_Container::Point_Container(size_type n_point, 
				 value_type c, 
				 value_type I,
				 value_type s)
     : last_idx_(n_point),
       ca_array_now_(n_point+2, c),
       ca_array_next_(n_point+2, c),
       ip3_array_now_(n_point+2, I),
       ip3_array_next_(n_point+2, I),
       s_array_now_(n_point+2, s),
       s_array_next_(n_point+2, s)
{
     assert(n_point != 0);

     /*
      * We actually store n_point+2 points :
      * - n_point actual points
      * - first point corresponds to uj-1*
      * - last point corresponds to ui+1*
      * (see Cell::compute_gap_junctions())
      */
}

// ====================

cellsim::size_type Point_Container::size() const
{
     // All arrays are guaranteed to have the same size
     return ca_array_now_.size();
}

void Point_Container::swap_times()
{
     std::swap(ca_array_now_, ca_array_next_);
     std::swap(ip3_array_now_, ip3_array_next_);
     std::swap(s_array_now_, s_array_next_);
}

// ====================
// Computing averages

namespace ba = boost::accumulators;
typedef ba::accumulator_set<double, ba::stats<ba::tag::mean> > mean_acc_type;

Point_Container::value_type Point_Container::c_avg() const
{
     mean_acc_type acc;
     std::for_each(ca_begin_now(), ca_end_now(),
		   [&acc] (const value_type& t)
		   {
			acc(t.value());
		   }
	  );

     return value_type::from_value(ba::mean(acc));
}

Point_Container::value_type Point_Container::I_avg() const
{
     mean_acc_type acc;
     std::for_each(ip3_begin_now(), ip3_end_now(),
		   [&acc] (const value_type& t)
		   {
			acc(t.value());
		   }
	  );

     return value_type::from_value(ba::mean(acc));
}

Point_Container::value_type Point_Container::s_avg() const
{
     mean_acc_type acc;
     std::for_each(s_begin_now(), s_end_now(),
		   [&acc] (const value_type& t)
		   {
			acc(t.value());
		   }
	  );

     return value_type::from_value(ba::mean(acc));
}

// ====================
// First/last element accessors

Point_Container::value_type Point_Container::ca_front_now() const
{
     return ca_array_now_[1];
}

Point_Container::value_type Point_Container::ca_back_now() const
{
     return ca_array_now_[last_idx_];
}


Point_Container::value_type Point_Container::ip3_front_now() const
{
     return ip3_array_now_[1];
}

Point_Container::value_type Point_Container::ip3_back_now() const
{
     return ip3_array_now_[last_idx_];
}


Point_Container::value_type Point_Container::s_front_now() const
{
     return s_array_now_[1];
}

Point_Container::value_type Point_Container::s_back_now() const
{
     return s_array_now_[last_idx_];
}


// ====================
// Boundaries accessors

void Point_Container::set_ca_bc(value_type front, value_type back)
{
     ca_array_now_.front() = front;
     ca_array_next_.front() = front;
     ca_array_now_.back() = back;
     ca_array_next_.back() = back;
}

void Point_Container::set_ip3_bc(value_type front, value_type back)
{
     ip3_array_now_.front() = front;
     ip3_array_next_.front() = front;
     ip3_array_now_.back() = back;
     ip3_array_next_.back() = back;
}

void Point_Container::set_s_bc(value_type front, value_type back)
{
     s_array_now_.front() = front;
     s_array_next_.front() = front;
     s_array_now_.back() = back;
     s_array_next_.back() = back;
}


Point_Container::value_type& Point_Container::ca_bc_front_now()
{
     return ca_array_now_.front();
}

Point_Container::value_type& Point_Container::ca_bc_back_now()
{
     return ca_array_now_.back();
}
Point_Container::value_type& Point_Container::ca_bc_front_next()
{
     return ca_array_next_.front();
}

Point_Container::value_type& Point_Container::ca_bc_back_next()
{
     return ca_array_next_.back();
}


Point_Container::value_type& Point_Container::ip3_bc_front_now()
{
     return ip3_array_now_.front();
}

Point_Container::value_type& Point_Container::ip3_bc_back_now()
{
     return ip3_array_now_.back();
}
Point_Container::value_type& Point_Container::ip3_bc_front_next()
{
     return ip3_array_next_.front();
}

Point_Container::value_type& Point_Container::ip3_bc_back_next()
{
     return ip3_array_next_.back();
}


Point_Container::value_type Point_Container::s_bc_front_now() const
{
     return s_array_now_.front();
}

Point_Container::value_type Point_Container::s_bc_back_now() const
{
     return s_array_now_.back();
}
Point_Container::value_type& Point_Container::s_bc_front_now()
{
     return s_array_now_.front();
}

Point_Container::value_type& Point_Container::s_bc_back_now()
{
     return s_array_now_.back();
}
Point_Container::value_type& Point_Container::s_bc_front_next()
{
     return s_array_next_.front();
}

Point_Container::value_type& Point_Container::s_bc_back_next()
{
     return s_array_next_.back();
}

// ====================
// Iterator accessors

// Ca
Point_Container::const_iterator Point_Container::ca_begin_now() const
{
     return ca_array_now_.begin() + 1;
}

Point_Container::const_iterator Point_Container::ca_end_now() const
{
     return ca_array_now_.end() - 1;
}

Point_Container::const_iterator Point_Container::ca_begin_next() const
{
     return ca_array_next_.begin() + 1;
}

Point_Container::const_iterator Point_Container::ca_end_next() const
{
     return ca_array_next_.end() - 1;
}

Point_Container::iterator Point_Container::ca_begin_next()
{
     return ca_array_next_.begin() + 1;
}

Point_Container::iterator Point_Container::ca_end_next()
{
     return ca_array_next_.end() - 1;
}

// IP3
Point_Container::const_iterator Point_Container::ip3_begin_now() const
{
     return ip3_array_now_.begin() + 1;
}

Point_Container::const_iterator Point_Container::ip3_end_now() const
{
     return ip3_array_now_.end() - 1;
}

Point_Container::const_iterator Point_Container::ip3_begin_next() const
{
     return ip3_array_next_.begin() + 1;
}

Point_Container::const_iterator Point_Container::ip3_end_next() const
{
     return ip3_array_next_.end() - 1;
}

Point_Container::iterator Point_Container::ip3_begin_next()
{
     return ip3_array_next_.begin() + 1;
}

Point_Container::iterator Point_Container::ip3_end_next()
{
     return ip3_array_next_.end() - 1;
}

// s
Point_Container::const_iterator Point_Container::s_begin_now() const
{
     return s_array_now_.begin() + 1;
}

Point_Container::const_iterator Point_Container::s_end_now() const
{
     return s_array_now_.end() - 1;
}

Point_Container::const_iterator Point_Container::s_begin_next() const
{
     return s_array_next_.begin() + 1;
}

Point_Container::const_iterator Point_Container::s_end_next() const
{
     return s_array_next_.end() - 1;
}

Point_Container::iterator Point_Container::s_begin_next()
{
     return s_array_next_.begin() + 1;
}

Point_Container::iterator Point_Container::s_end_next()
{
     return s_array_next_.end() - 1;
}
