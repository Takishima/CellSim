/*
 * (c) All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
 * Switzerland, Laboratory of Cell Biophysics, Nguyen Damien, 2013
 * See the LICENSE file for more details.
 */
 
#ifndef CELL_IMPL_HPP_INCLUDED
#define CELL_IMPL_HPP_INCLUDED

/*
 * This is not pretty, but makes the code much more readable.
 *
 * What this macro basically does is to remove the ambiguity as to
 * which method of Cell to call (since there is no other way I know of...)
 */
/*!
 * Return a pointer to the non-const version of a calcium accessor of a Cell
 * \param x Address of a member function
 * \return Pointer to member function
 */
#define GET_METHOD_PTR_C(x)					\
     std::mem_fn(static_cast<concentration_ut&(Cell::*)()>(x))
//! Overload of \ref GET_METHOD_PTR_C for I
#define GET_METHOD_PTR_I(x) GET_METHOD_PTR_C(x)
//! Overload of \ref GET_METHOD_PTR_C for s
#define GET_METHOD_PTR_S(x) GET_METHOD_PTR_C(x)
//! Overload of \ref GET_METHOD_PTR_C for v
#define GET_METHOD_PTR_V(x)						\
     std::mem_fn(static_cast<electric_potential_ut&(Cell::*)()>(x))
//! Overload of \ref GET_METHOD_PTR_C for w
#define GET_METHOD_PTR_W(x) \
     std::mem_fn(static_cast<double&(Cell::*)()>(x))

/*!
 * Synctatic sugar to transform the model equations of motion to the
 * appropriate functor for the integrator methods.
 * \param name Name of function to bind
 * \return Functor
 */
#define GET_EQN_PTR(name) std::bind(name, _1, std::cref(model_constants), _2)

#endif //CELL_IMPL_HPP_INCLUDED
