#ifndef MPM_MPM_EXPLICIT_USL_H_
#define MPM_MPM_EXPLICIT_USL_H_

#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "container.h"
#include "mpm.h"
#include "mpm_explicit.h"
#include "particle.h"

namespace mpm {

//! MPMExplicitUSl class
//! \brief Explicit one phase mpm with USL
//! \details A single-phase explicit MPM with Update Stress Last
//! \tparam Tdim Dimension
template <unsigned Tdim>
class MPMExplicitUSL : public MPMExplicit<Tdim> {
 public:
  //! Constructor
  MPMExplicitUSL(std::unique_ptr<IO>&& io);

  //! Solve
  bool solve() override;

 protected:
  // Generate a unique id for the analysis
  using mpm::MPMExplicit<Tdim>::uuid_;
  //! Time step size
  using mpm::MPMExplicit<Tdim>::dt_;
  //! Number of steps
  using mpm::MPMExplicit<Tdim>::nsteps_;
  //! Output steps
  using mpm::MPMExplicit<Tdim>::output_steps_;
  //! A unique ptr to IO object
  using mpm::MPMExplicit<Tdim>::io_;
  //! JSON analysis object
  using mpm::MPMExplicit<Tdim>::analysis_;
  //! JSON post-process object
  using mpm::MPMExplicit<Tdim>::post_process_;
  //! Logger
  using mpm::MPMExplicit<Tdim>::console_;

  //! Gravity
  using mpm::MPMExplicit<Tdim>::gravity_;
  //! Mesh object
  using mpm::MPMExplicit<Tdim>::meshes_;
  //! Materials
  using mpm::MPMExplicit<Tdim>::materials_;

};  // MPMExplicitUSl class
}  // namespace mpm

#include "mpm_explicit_usl.tcc"

#endif  // MPM_MPM_EXPLICIT_USL_H_