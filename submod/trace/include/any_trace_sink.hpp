#pragma once

#include <pemu/trace/trace.hpp>

#include <concepts>
#include <memory>
#include <type_traits>
#include <utility>

namespace pemu::trace {

/**
 * @brief Owning runtime trace sink for non-template public APIs.
 *
 * Concrete sinks should normally be used through the TraceSink concept so
 * calls can be inlined. AnyTraceSink is reserved for boundaries such as the
 * non-template Simulation facade, where the concrete sink type cannot be part
 * of the class type. Copies share one underlying sink and its buffered state.
 */
class AnyTraceSink {
 public:
  AnyTraceSink() noexcept = default;

  template <typename Sink>
    requires(!std::same_as<std::remove_cvref_t<Sink>, AnyTraceSink>) &&
            TraceSink<std::remove_cvref_t<Sink>> &&
            std::constructible_from<std::remove_cvref_t<Sink>, Sink>
  AnyTraceSink(Sink&& sink)
      : sink_(std::make_shared<Model<std::remove_cvref_t<Sink>>>(
            std::forward<Sink>(sink))) {}

  void operator()(const TraceEvent& event) const noexcept {
    if (sink_) {
      sink_->emit(event);
    }
  }

  void flush() const noexcept {
    if (sink_) {
      sink_->flush();
    }
  }

  [[nodiscard]] explicit operator bool() const noexcept {
    return static_cast<bool>(sink_);
  }

 private:
  struct Interface {
    virtual ~Interface() = default;
    virtual void emit(const TraceEvent& event) noexcept = 0;
    virtual void flush() noexcept = 0;
  };

  template <TraceSink Sink>
  struct Model final : Interface {
    template <typename Value>
      requires std::constructible_from<Sink, Value>
    explicit Model(Value&& value) : sink(std::forward<Value>(value)) {}

    void emit(const TraceEvent& event) noexcept override { sink(event); }
    void flush() noexcept override { sink.flush(); }

    Sink sink;
  };

  std::shared_ptr<Interface> sink_;
};

static_assert(TraceSink<AnyTraceSink>);

}  // namespace pemu::trace
