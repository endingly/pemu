#pragma once

#include <pemu/output/dump/types.hpp>

#include <memory>

namespace pemu::output::dump {

/** @brief Owns one ordered visualization-dump time series. */
class ISeries {
 public:
  virtual ~ISeries() = default;

  [[nodiscard]] virtual OutputRecord append(const Request& request) = 0;
  [[nodiscard]] virtual OutputRecord finish() = 0;
};

/** @brief Factory for independent visualization-dump series. */
class IWriter {
 public:
  virtual ~IWriter() = default;

  [[nodiscard]] virtual std::unique_ptr<ISeries> openSeries(
      const SeriesRequest& request) const = 0;

  [[nodiscard]] OutputRecord write(const Request& request) const {
    auto series = openSeries({.mesh = request.mesh,
                              .path = request.path,
                              .overwrite = request.overwrite});
    (void)series->append(request);
    return series->finish();
  }
};

}  // namespace pemu::output::dump
