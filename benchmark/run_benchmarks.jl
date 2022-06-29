# To run it locally, BenchmarkCI should be added to root project
using BenchmarkCI
on_CI = haskey(ENV, "GITHUB_ACTIONS")

judgement = BenchmarkCI.judge()
if on_CI
    BenchmarkCI.postjudge()
else
    ciresult = BenchmarkCI.CIResult(; judgement=judgement)
    output_file = "benchmark_result.md"
    @info "benchmark finished" output_file
    open(output_file, "w") do io
        BenchmarkCI.printcommentmd(io, ciresult)
    end
end
