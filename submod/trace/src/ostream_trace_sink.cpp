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

    if (stream_->good() && shouldFlushAfter(event)) {
      stream_->flush();
    }
    failed_ = failed_ || !stream_->good();
  } catch (...) {
    failed_ = true;
  }
}

bool OstreamTraceSink::shouldFlushAfter(const TraceEvent& event) noexcept {
  if (event.domain != DiagDomain::simulation) {
    return false;
  }

  if (event.name == "step.completed") {
    ++completed_steps_since_flush_;
    if (completed_steps_since_flush_ == 2) {
      completed_steps_since_flush_ = 0;
      return true;
    }
    return false;
  }

  if (event.name == "run.completed" || event.name == "run.failed") {
    completed_steps_since_flush_ = 0;
    return true;
  }

  return false;
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
