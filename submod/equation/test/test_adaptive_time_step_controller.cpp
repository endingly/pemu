#include <gtest/gtest.h>

#include <pemu/equation/time_integration/adaptive_time_step_controller.hpp>

#include <stdexcept>

namespace pemu::equation::time_integration::test {

TEST(AdaptiveTimeStepControllerTest, AppliesSafetyAndGrowthLimits) {
  AdaptiveTimeStepController controller(
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 0.5, .max_growth = 1.5});

  const auto proposal = controller.propose(1.0, 0.5, 0.1, 1.0);

  EXPECT_NEAR(proposal.transport_limit, 1.0, 1e-15);
  EXPECT_NEAR(proposal.positivity_limit, 0.5, 1e-15);
  EXPECT_NEAR(proposal.stability_limit, 0.5, 1e-15);
  EXPECT_NEAR(proposal.dt, 0.15, 1e-15);
}

TEST(AdaptiveTimeStepControllerTest, TruncatesFinalStepToRemainingTime) {
  AdaptiveTimeStepController controller(
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 0.5, .max_growth = 2.0});

  const auto proposal = controller.propose(1.0, 1.0, 0.0, 0.03);

  EXPECT_NEAR(proposal.dt, 0.03, 1e-15);
}

TEST(AdaptiveTimeStepControllerTest, RejectsNoPositiveAdmissibleTimeStep) {
  AdaptiveTimeStepController controller(
      {.safety = 0.9, .min_dt = 1e-8, .max_dt = 0.5, .max_growth = 2.0});

  EXPECT_THROW((void)controller.propose(1.0, 0.0, 0.0, 1.0),
               std::runtime_error);
}

}  // namespace pemu::equation::time_integration::test
