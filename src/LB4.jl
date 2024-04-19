module LB4

struct LoopInfo
  times::Vector{Float64}
  nb_threads::UInt64
end

function Base.show(io::IO, li::LoopInfo)
  digits = 3
  mean_time = sum(li.times) / li.nb_threads
  sd_time = sqrt(sum(map(x -> x^2, li.times .- mean_time)) / (li.nb_threads - 1))
  cov = round(sd_time / mean_time, digits=digits)
  max_time = maximum(li.times) 
  bound = round(1.96 * sd_time / sqrt(li.nb_threads), digits=digits)
  percentage_imbalance = round((li.nb_threads / (li.nb_threads-1)) * ((max_time - mean_time) / max_time) * 100, digits=digits)
  print(io, "$(round(mean_time, digits=digits))s ± $(bound) (95% CI), c.o.v.: $(cov), p.i.: $(percentage_imbalance)%")
end

function threading_run(fun)
    ccall(:jl_enter_threaded_region, Cvoid, ())
    n = Threads.threadpoolsize()
    tid_offset = Threads.threadpoolsize(:interactive)
    tasks = Vector{Task}(undef, n)
    times = Vector{Float64}(undef, n)
    start_time = time_ns()
    for i = 1:n
        t = Task(() -> fun()) # pass in tid
        t.sticky = true
        ccall(:jl_set_task_tid, Cint, (Any, Cint), t, tid_offset + i-1)
        tasks[i] = t
        schedule(t)
    end
    for i = 1:n
        Base._wait(tasks[i])
        times[i] = (fetch(tasks[i]) - start_time) * 1e-9
    end
    ccall(:jl_exit_threaded_region, Cvoid, ())
    failed_tasks = filter!(istaskfailed, tasks)
    if !isempty(failed_tasks)
        throw(CompositeException(map(TaskFailedException, failed_tasks)))
    end
    LoopInfo(times, n)
end

function get_chunk(sched::Symbol, n::Int64, p::Int64)
  if sched == :static
    (i::Int64) -> ceil(Int64, n / p)
  elseif sched == :dynamic
    (i::Int64) -> 1
  elseif sched == :gss
    (i::Int64) -> ceil(Int64, (n / p) * ((p - 1) / p)^i)
  elseif sched == :fac2
    (i::Int64) -> ceil(Int64, ((n / p) * 0.5^(floor(Int64, i/p) + 1)))
  elseif sched == :tss
    k_0 = ceil(Int64, n / (2 * p))
    k_s_1 = 1
    s = ceil(Int64, 2 * n / (k_0 + k_s_1))
    k = floor(Int64, (k_0 - k_s_1) / (s - 1))
    (i::Int64) -> (k_0 - i * k)
  elseif sched == :fiss
    b = ceil(Int64, log2(n))
    k_0 = ceil(Int64, n / ((2 + b) * p))
    k = ceil(Int64, (2 * n * (1 - (b / (2 + b)))) / (p * b * (b - 1)))
    (i::Int64) -> (k_0 + i * k) 
  elseif sched == :viss
    b = ceil(Int64, log2(n))
    k_0 = ceil(Int64, n / ((2 + b) * p))
    (i::Int64) -> k_0 * ceil(Int64, (1 - (1/2)^(i % p)) / (1/2))
  elseif sched == :rnd
    k_min = ceil(Int64, n / (100 * p))
    k_max = ceil(Int64, n / (2 * p))
    (i::Int64) -> rand((k_min:k_max))
  elseif sched == :auto
    phi = (1 + sqrt(5)) / 2
    f = floor(Int64, log2(n / p) / phi)
    k = floor(Int64, n / (2 * p * (2^f)))
    (i::Int64) -> k
  else
    println("Unknown sched: ", sched)
    throw("Unknown `sched`")
  end
end

struct LogInfo
  start_ts::UInt64
  end_ts::UInt64
  thread_id::Int
  start_iter::Int64
  end_iter::Int64
end

Base.show(io::IO, li::LogInfo) = print(io, "$(li.thread_id), $(li.start_iter), $(li.end_iter), $(li.start_ts), $(li.end_ts)")

get_str(li::LogInfo, sched::Symbol) = "$(li.thread_id), $(li.start_iter), $(li.end_iter), $(li.start_ts), $(li.end_ts), $(sched)\n"

macro lbthreads_log(args...)
    log_prefix = "lb4jl_thread_" 
    if length(args) == 3
      sched, log_prefix, loop = args
    else
      sched, loop = args
    end
    iterator = loop.args[1]
    body = loop.args[2]
    index_variable = iterator.args[1]
    range = iterator.args[2]
    quote
      local thread_f
      let range = $(esc(range))
        lenr = length(range)
        let sched_v = $(esc(sched))
        let log_prefix_v = $(esc(log_prefix))

        p = Threads.threadpoolsize()
        chunk_f = typeof(sched_v) == Symbol ? get_chunk(sched_v, lenr, p) : eval(:($sched_v($lenr, $p)))

        queue_index = Threads.Atomic{Int}(0)
        fi = firstindex(range)
        start = Threads.Atomic{Int}(fi)

        function thread_f()
          tasks_log = LogInfo[]
          while true
            index = Threads.atomic_add!(queue_index, 1)
            block = chunk_f(index)
            start_iter = Threads.atomic_add!(start, block + 1)
            end_iter = start_iter + block

            if start_iter >= fi + lenr ||  end_iter < start_iter
              #println("bogus: ($(start_iter), $(end_iter))")
              break
            end
            if end_iter >= fi + lenr
              end_iter = lenr
            end
            start_ts = time_ns()
            for i in start_iter:end_iter
              local $(esc(index_variable)) = @inbounds i
              $(esc(body))
            end
            end_ts = time_ns()
            push!(tasks_log, LogInfo(start_ts, end_ts, Threads.threadid(), start_iter, end_iter))
          end
          outfile_name = "$(log_prefix_v)$(Threads.threadid()).csv" 
          open(outfile_name, "w") do file
            for task in tasks_log
              write(file, get_str(task, sched_v))
            end
          end
          time_ns()
        end
      end
      end
      end
      threading_run(thread_f)
    end
end

macro lbthreads(args...)
    sched, loop = args
    iterator = loop.args[1]
    body = loop.args[2]
    index_variable = iterator.args[1]
    range = iterator.args[2]
    quote
      local thread_f
      let range = $(esc(range))
        lenr = length(range)
        let sched_v = $(esc(sched))

        p = Threads.threadpoolsize()
        chunk_f = typeof(sched_v) == Symbol ? get_chunk(sched_v, lenr, p) : eval(:($sched_v($lenr, $p)))

        queue_index = Threads.Atomic{Int}(0)
        fi = firstindex(range)
        start = Threads.Atomic{Int}(fi)

        function thread_f()
          while true
            index = Threads.atomic_add!(queue_index, 1)
            block = chunk_f(index)
            start_iter = Threads.atomic_add!(start, block + 1)
            end_iter = start_iter + block
            if start_iter >= fi + lenr ||  end_iter < start_iter
              #println("bogus")
              break
            end
            if end_iter >= fi + lenr
              end_iter = lenr
            end
            for i in start_iter:end_iter
              local $(esc(index_variable)) = @inbounds i
              $(esc(body))
            end
          end
          time_ns()
        end
      end
    end
      threading_run(thread_f)
    end
end

end
