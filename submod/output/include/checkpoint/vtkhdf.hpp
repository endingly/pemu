#pragma once

#include <pemu/output/checkpoint/i_reader.hpp>
#include <pemu/output/checkpoint/i_writer.hpp>

namespace pemu::output::checkpoint {

inline constexpr std::uint32_t vtkhdf_format_version = 2;

/** @brief Writes exact restart fields into a versioned VTKHDF checkpoint. */
class VtkHdfWriter final : public IWriter {
 public:
  [[nodiscard]] OutputRecord write(const WriteRequest& request) const override;
};

/** @brief Inspects and restores versioned VTKHDF checkpoints. */
class VtkHdfReader final : public IReader {
 public:
  [[nodiscard]] Manifest inspect(
      const std::filesystem::path& path) const override;

  [[nodiscard]] OutputRecord restore(
      const std::filesystem::path& path,
      const RestoreRequest& request) const override;
};

}  // namespace pemu::output::checkpoint
