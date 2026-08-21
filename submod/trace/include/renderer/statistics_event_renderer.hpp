#pragma once

#include <pemu/trace/trace.hpp>

#include <cstdint>
#include <iosfwd>
#include <optional>
#include <string>
#include <vector>

namespace pemu::trace::renderer {

class StatisticsEventRenderer {
 public:
  void render(std::ostream& stream, const TraceEvent& event);
  void flush(std::ostream& stream);

 private:
  struct Row {
    std::string field;
    std::string unit;
    std::string minimum;
    std::string maximum;
    std::string mean;
    std::string rms;
    std::string integral;
    Severity severity{Severity::Debug};
    EventKind kind{EventKind::Trace};
    std::uint64_t non_finite{};
    std::uint64_t negative{};
    std::uint64_t invalid_weight{};
  };

  struct Group {
    bool active{false};
    std::optional<std::uint64_t> step;
    std::string time;
    std::optional<std::uint64_t> cell_samples;
    std::optional<std::uint64_t> face_samples;
    std::string volume_semantics;
    std::string physical_volume;
    std::string physical_face_area;
    std::vector<Row> rows;
  };

  void renderGroup(std::ostream& stream);

  Group group_;
  bool has_output_{false};
  bool metadata_written_{false};
};

}  // namespace pemu::trace::renderer
