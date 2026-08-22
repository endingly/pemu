#include <pemu/trace/ostream_trace_sink.hpp>

#include <ostream>

namespace pemu::trace {

OstreamTraceSink::OstreamTraceSink(std::ostream& stream) noexcept
    : stream_(&stream) {}

void OstreamTraceSink::operator()(const TraceEvent& event) noexcept {
  try {
    if (event.output_channel == OutputChannel::Statistics) {
      statistics_renderer_.render(*stream_, event);
      failed_ = failed_ || !stream_->good();
      return;
    }

    statistics_renderer_.flush(*stream_);
    const bool write_header = !event_header_written_;
    if (event.kind == EventKind::Diagnostic) {
      diagnostic_renderer_.render(*stream_, event, write_header);
    } else {
      ordinary_renderer_.render(*stream_, event, write_header);
    }
    event_header_written_ = true;
    failed_ = failed_ || !stream_->good();
  } catch (...) {
    failed_ = true;
  }
}

void OstreamTraceSink::flush() noexcept {
  try {
    statistics_renderer_.flush(*stream_);
    stream_->flush();
    failed_ = failed_ || !stream_->good();
  } catch (...) {
    failed_ = true;
  }
}

bool OstreamTraceSink::failed() const noexcept {
  return failed_;
}

}  // namespace pemu::trace
