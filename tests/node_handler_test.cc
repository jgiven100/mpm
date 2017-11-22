// No-debase test
#include <limits>

#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "node.h"
#include "node_handler.h"

//! \brief Check node handler class for 2D case
TEST_CASE("Node handler is checked for 2D case", "[node][2D]") {
  const unsigned Dim = 2;

  // Node 1
  mpm::Index id1 = 0;
  Eigen::Vector2d coords;
  coords.setZero();
  auto node1 = std::make_shared<mpm::Node<Dim>>(id1, coords);

  // Node 2
  mpm::Index id2 = 1;
  auto node2 = std::make_shared<mpm::Node<Dim>>(id2, coords);

  // Node handler
  auto nodehandler = std::make_shared<mpm::NodeHandler<Dim>>();

  // Check insert node
  SECTION("Check insert node functionality") {
    // Insert node 1 and check status
    bool status1 = nodehandler->insert_node(node1);
    REQUIRE(status1 == true);
    // Insert node 2 and check status
    bool status2 = nodehandler->insert_node(node2);
    REQUIRE(status2 == true);
    // Check size of node hanlder
    REQUIRE(nodehandler->nnodes() == 2);
  }

  // Check iterator
  SECTION("Check begin node iterator") {
    // Insert node 1
    nodehandler->insert_node(node1);
    // Insert node 2
    nodehandler->insert_node(node2);
    // Check size of node hanlder
    std::size_t counter = 0;
    for (auto itr = nodehandler->nodes_begin(); itr != nodehandler->nodes_end();
         ++itr) {
      auto coords = ((*itr).second)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i) REQUIRE(coords[i] == 0);
      ++counter;
    }
    // Iterate over nodes and check if the number of nodes is good
    REQUIRE(counter == 2);
  }
}

//! \brief Check node handler class for 3D case
TEST_CASE("Node handler is checked for 3D case", "[nodehandler][3D]") {
  const unsigned Dim = 3;

  // Node 1
  mpm::Index id1 = 0;
  Eigen::Vector3d coords;
  coords.setZero();
  auto node1 = std::make_shared<mpm::Node<Dim>>(id1, coords);

  // Node 2
  mpm::Index id2 = 1;
  auto node2 = std::make_shared<mpm::Node<Dim>>(id2, coords);

  // Node handler
  auto nodehandler = std::make_shared<mpm::NodeHandler<Dim>>();

  // Check insert node
  SECTION("Check insert node functionality") {
    // Insert node 1 and check status
    bool status1 = nodehandler->insert_node(node1);
    REQUIRE(status1 == true);
    // Insert node 2 and check status
    bool status2 = nodehandler->insert_node(node2);
    REQUIRE(status2 == true);
    // Check size of node hanlder
    REQUIRE(nodehandler->nnodes() == 2);
  }

  // Check iterator
  SECTION("Check begin node iterator") {
    // Insert node 1
    nodehandler->insert_node(node1);
    // Insert node 2
    nodehandler->insert_node(node2);
    // Check size of node hanlder
    std::size_t counter = 0;
    for (auto itr = nodehandler->nodes_begin(); itr != nodehandler->nodes_end();
         ++itr) {
      auto coords = ((*itr).second)->coordinates();
      // Check if coordinates for each node is zero
      for (unsigned i = 0; i < coords.size(); ++i) REQUIRE(coords[i] == 0);
      ++counter;
    }
    // Iterate over nodes and check if the number of nodes is good
    REQUIRE(counter == 2);
  }
}
