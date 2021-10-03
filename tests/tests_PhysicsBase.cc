#include <gtest/gtest.h>

#include "PhysicsBase.hh"

#include <cmath>

using namespace CC6::literals;

TEST(tests_PhysicsBase, literals) {
  auto a1 = 0.0_rad;
  auto a2 = 90.0_rad;
  auto a3 = 180.0_rad;
  auto a4 = -90.0_rad;
  auto a5 = -360.0_rad;

  ASSERT_DOUBLE_EQ(a1, 0.);
  ASSERT_DOUBLE_EQ(a2, M_PI / 2);
  ASSERT_DOUBLE_EQ(a3, M_PI);
  ASSERT_DOUBLE_EQ(a4, -M_PI / 2);
  ASSERT_DOUBLE_EQ(a5, -M_PI * 2);

  auto a11 = 0_rad;
  auto a12 = 90_rad;
  auto a13 = 180_rad;
  auto a14 = -90_rad;
  auto a15 = -360_rad;

  ASSERT_DOUBLE_EQ(a11, 0.);
  ASSERT_DOUBLE_EQ(a12, M_PI / 2);
  ASSERT_DOUBLE_EQ(a13, M_PI);
  ASSERT_DOUBLE_EQ(a14, -M_PI / 2);
  ASSERT_DOUBLE_EQ(a15, -M_PI * 2);
}

TEST(tests_PhysicsBase, cs_gamma_E) {
  auto E = 100.0;

  auto e1 = CC6::ComptonScatteringGammaE(0_rad, 100);
  ASSERT_EQ(e1, E);

  auto e2 = CC6::ComptonScatteringGammaE(360_rad, 100);
  ASSERT_EQ(e2, E);

  auto e3 = CC6::ComptonScatteringGammaE(0.0_rad, 100);
  ASSERT_EQ(e3, E);

  auto e4 = CC6::ComptonScatteringGammaE(360.0_rad, 100);
  ASSERT_EQ(e4, E);
}
