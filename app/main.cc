#include <fmt/core.h>
#include <mdspan/mdspan.hpp>
#include <taskflow/taskflow.hpp>

#include <vector>

constexpr auto N = 100;

int main() {
  tf::Executor exec;
  tf::Taskflow taskflow;

  std::vector<double> a(N * N);

  Kokkos::mdspan mda(a.data(), N, N);

  auto [A, B, C, D] = taskflow.emplace(
      []() { fmt::print("task A\n"); }, []() { fmt::print("task B\n"); },
      []() { fmt::print("task C\n"); }, []() { fmt::print("task D\n"); });

  A.precede(B, C);
  D.succeed(B, C);

  exec.run(taskflow).wait();

  return 0;
}
