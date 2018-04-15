#include "koenigsberger_activator_base.hpp"

namespace activator = cellsim::model_koenigsberger::activator;

activator::Activator_Base::Activator_Base(size_type start_column,
				size_type end_column)
     : start_column_(start_column), end_column_(end_column)
{}

bool activator::Activator_Base::do_activation(size_type index) const
{
     return (start_column_ <= index && index < end_column_);
}

