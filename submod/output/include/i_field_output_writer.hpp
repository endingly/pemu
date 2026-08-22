#pragma once

#include <pemu/output/output_types.hpp>

#include <memory>

namespace pemu::output {

/** @brief Owns the state of one ordered field-output time series. */
class IFieldOutputSeries {
 public:
  /** @brief Destroys an output-series session through its interface. */
  virtual ~IFieldOutputSeries() = default;

  /**
   * @brief Captures the next synchronized field state in strictly increasing
   * time order.
   * @param request Snapshot fields and temporal stamp.
   * @return Logical record identifying the series file and captured state.
   */
  [[nodiscard]]
  virtual OutputRecord append(const FieldDumpRequest& request) = 0;

  /**
   * @brief Finalizes the series and makes its output file readable.
   * @return Record identifying the completed file and its last state.
   */
  [[nodiscard]]
  virtual OutputRecord finish() = 0;
};

/** @brief Factory for independent field-output series. */
class IFieldOutputWriter {
 public:
  virtual ~IFieldOutputWriter() = default;

  /**
   * @brief Opens an ordered field-output series.
   * @param request Mesh, destination path, and overwrite policy.
   * @return A stateful series that must be finished after its last snapshot.
   */
  [[nodiscard]]
  virtual std::unique_ptr<IFieldOutputSeries> openSeries(
      const FieldSeriesRequest& request) const = 0;

  /**
   * @brief Writes one standalone snapshot through a one-state series.
   * @param request Complete standalone field-dump request.
   * @return Record identifying the completed output file.
   */
  [[nodiscard]]
  OutputRecord write(const FieldDumpRequest& request) const {
    auto series = openSeries({.mesh = request.mesh,
                              .path = request.path,
                              .overwrite = request.overwrite});
    (void)series->append(request);
    return series->finish();
  }
};

}  // namespace pemu::output
