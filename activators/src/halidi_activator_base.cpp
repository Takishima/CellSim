#include "halidi_activator_base.hpp"

namespace activator = cellsim::model_halidi::activator;

activator::Activator_Base::Activator_Base(size_type cell_idx)
     : cell_idx_(cell_idx)
{}

bool activator::Activator_Base::do_activation(size_type i) const
{
     return i == cell_idx_;
}
