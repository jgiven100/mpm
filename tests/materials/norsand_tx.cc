#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include <cmath>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "materials/material.h"
#include "node.h"
#include "particle.h"

//! Check NorSand class in 3D
TEST_CASE("NorSand is checked in 3D", "[material][NorSand][3D]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  const unsigned Dim = 3;

  // Add particle
  mpm::Index pid = 0;
  Eigen::Matrix<double, Dim, 1> coords;
  coords.setZero();
  auto particle = std::make_shared<mpm::Particle<Dim>>(pid, coords);

  // Initialise material
  Json jmaterial;

  // Testing
  jmaterial["density"] = 2000.;
  jmaterial["poisson_ratio"] = 0.15;
  jmaterial["reference_pressure"] = 1.0E+5;
  jmaterial["friction_cs"] = 31.125;
  jmaterial["N"] = 0.35;
  jmaterial["lambda"] = 0.04;
  jmaterial["kappa"] = 0.00589;
  jmaterial["gamma"] = 0.691;
  jmaterial["chi"] = 4.;
  jmaterial["hardening_modulus"] = 100.0;
  jmaterial["void_ratio_initial"] = 0.613067;
  jmaterial["p_image_initial"] = 73575.9;
  jmaterial["bond_model"] = false;
  jmaterial["p_cohesion_initial"] = 2.0E+3;
  jmaterial["p_dilation_initial"] = 5.0E+3;
  jmaterial["m_cohesion"] = 2.0;
  jmaterial["m_dilation"] = 2.0;
  jmaterial["m_modulus"] = 1000;

  double pi_constant = M_PI;

  // Check triaxial drain test
  SECTION("NorSand check drained stresses") {
    unsigned id = 0;
    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "NorSand3D", std::move(id), jmaterial);
    REQUIRE(material->id() == 0);


    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    stress(0) = -200000;
    stress(1) = -200000;
    stress(2) = -200000-1;


    // Check dmatrix
    const double poisson_ratio_ = jmaterial["poisson_ratio"];
    const double kappa = jmaterial["kappa"];
    const double e_init = jmaterial["void_ratio_initial"];
    const double p_dilation = jmaterial["m_dilation"];
    const double p_cohesion = jmaterial["m_cohesion"];
    const double m_modulus = jmaterial["m_modulus"];
    const double p = (stress(0) + stress(1) + stress(2)) / 3.0;
    const double bulk_modulus_ =
        (1. + e_init) / kappa * p + m_modulus * (p_cohesion + p_dilation);
    const double G = 3. * bulk_modulus_ * (1. - 2. * poisson_ratio_) /
                     (2.0 * (1. + poisson_ratio_));
    const double a1 = bulk_modulus_ + (4.0 / 3.0) * G;
    const double a2 = bulk_modulus_ - (2.0 / 3.0) * G;


    // Initialise dstrain (UNDRAINED)
    unsigned multiplier = 50;
    mpm::Material<Dim>::Vector6d dstrain;
    dstrain.setZero();
    dstrain(2) = -0.0005 / multiplier;
    dstrain(0) = -0.5 * dstrain(2);
    dstrain(1) = -0.5 * dstrain(2);


    // Compute updated stress
    mpm::dense_map state_vars = material->initialise_state_variables();
    std::ofstream myfile;
    myfile.open("norsand_tx_undrained.txt");
    double axial_strain = 0.;


    // Write initial states
    myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
           << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
           << axial_strain << '\t' << (state_vars).at("void_ratio") << '\n';

    // Loop
    for (unsigned i = 0; i < multiplier * 1000 - 1; ++i) {
      axial_strain += dstrain(2);
      stress = material->compute_stress(stress, dstrain, particle.get(),
                                        &state_vars);

      if (i % multiplier == 0) {
         myfile << stress(0) << '\t' << stress(1) << '\t' << stress(2) << '\t'
                << stress(3) << '\t' << stress(4) << '\t' << stress(5) << '\t'
                << axial_strain << '\t' << (state_vars).at("void_ratio") <<'\n';
      }
    }
    myfile.close();
  }
}