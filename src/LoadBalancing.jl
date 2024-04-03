module LoadBalancing

function threading_run(fun, static)
    ccall(:jl_enter_threaded_region, Cvoid, ())
    n = Threads.threadpoolsize()
    tid_offset = Threads.threadpoolsize(:interactive)
    tasks = Vector{Task}(undef, n)
    for i = 1:n
        t = Task(() -> fun()) # pass in tid
        t.sticky = static
        if static
            ccall(:jl_set_task_tid, Cint, (Any, Cint), t, tid_offset + i-1)
        else
            # TODO: this should be the current pool (except interactive) if there
            # are ever more than two pools.
            @assert ccall(:jl_set_task_threadpoolid, Cint, (Any, Int8), t, _sym_to_tpid(:default)) == 1
        end
        tasks[i] = t
        schedule(t)
    end
    for i = 1:n
        Base._wait(tasks[i])
    end
    ccall(:jl_exit_threaded_region, Cvoid, ())
    failed_tasks = filter!(istaskfailed, tasks)
    if !isempty(failed_tasks)
        throw(CompositeException(map(TaskFailedException, failed_tasks)))
    end
end

function get_chunk_dynamic(j::Int64, first_index::Int64, max_len::Int64)::Tuple{Int, Int}
  (first_index + (j - 1), first_index + (j - 1))
end

function get_chunk_static(j::Int64, first_index::Int64, max_len::Int64, block_size::Int64)::Tuple{Int, Int}
  (first_index + (j - 1) * block_size, first_index + min(max_len, j * block_size - 1))
end

function get_chunk_gss(j::Int64, first_index::Int64, max_len::Int64, p::Int64)::Tuple{Int, Int}
  offset = max_len - max_len * ((p-1)/p)^(j-1)
  k = (max_len / p) * ((p-1)/p)^(j-1)
  (floor(first_index + offset), floor(first_index + offset + k - 1))
end

function get_chunk_tss(j::Int64, first_index::Int64, max_len::Int64, f::Int64, l::Int64, S::Int64, delta::Int64)::Tuple{Int, Int}
  i_start = first_index + floor(((j - 1) * (2 * f - (j - 2) * delta)) / 2)# + 1
  i_end   = first_index + floor((j * (2 * f - (j - 1) * delta)) / 2)

  (i_start, (i_end - 1 <= max_len) ? i_end - 1 : max_len)
end

function get_chunk_fac2(j::Int64, first_index::Int64, max_len::Int64, p::Int64)::Tuple{Int, Int}
  round = (j - 1) ÷ p
  block_size = ceil(Int64, (0.5)^(round + 1) * max_len / p)
  start = first_index + sum(p * (ceil(Int64, (0.5)^(i + 1) * max_len / p)) for i in 0:(round - 1); init=0) + ((j - 1) % p) * block_size

  #println("Thread id: ", Threads.threadid(), ", round: ", round, ", block size: ", block_size, ", start: ", start)

  (start, start + block_size - 1)
end

function get_chunk_fiss(j::Int64, first_index::Int64, block_size::Int64)::Tuple{Int, Int}
  start = first_index + (j - 1) * block_size 
  (start, start + block_size - 1)
end

# function get_chunk_viss(j::Int64, first_index::Int64, block_size::Int64, p::Int64)::Tuple{Int, Int}
#   round = (j - 1) ÷ p
# 
#   start = first_index + (j - 1) * block_size 
#   (start, start + block_size - 1)
# end



function get_chunk(first_index::Int64, max_len::Int64, sched::Symbol)
  p = Threads.threadpoolsize()
  if sched == :static
    ((i::Int64) -> get_chunk_static(i, first_index, max_len, div(max_len, p)), p)
  elseif sched == :dynamic
    ((i::Int64) -> get_chunk_dynamic(i, first_index, max_len), max_len)
  elseif sched == :gss
    ((i::Int64) -> get_chunk_gss(i, first_index, max_len, p), ceil(log(p / max_len) / log((p-1)/p) + 1))
  elseif sched == :tss
    f = ceil(Int64, max_len / (2 * Threads.threadpoolsize()))
    l = 1
    S = ceil(Int64, 2 * max_len / (f + l))
    delta = floor(Int64, (f - l) / (S - 1))
    ((i::Int64) -> get_chunk_tss(i, first_index, max_len, f, l, S, delta), S)
  elseif sched == :fac2
    ((i::Int64) -> get_chunk_fac2(i, first_index, max_len, p), max_len)
  elseif sched == :fiss
    B = 1
    block_size = ceil((2 * max_len * (1 - (B / (2 + B)))) / (p * B * (B - 1)))
    ((i::Int64) -> get_chunk_fiss(i, first_index, block_size), max_len) # TODO: max_len
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



macro lbthreads(args...)
    sched, loop = args
    iterator = loop.args[1]
    body = loop.args[2]
    index_variable = iterator.args[1]
    range = iterator.args[2]
    quote
      local thread_f
      let range = $(esc(range))
        let sched_v = $(esc(sched))
        lenr = length(range)
        fi = firstindex(range)
        chunk_f, max_chunks = get_chunk(fi, lenr, sched_v)
        queue_index = Threads.Atomic{Int}(1)

        function thread_f()
          tasks_log = LogInfo[]
          while true
            index = Threads.atomic_add!(queue_index, 1)
            if index > max_chunks
              break
            end
            start_iter, end_iter = chunk_f(index)
            #end_iter = min(end_iter, fi + lenr)
            if start_iter >= fi + lenr || end_iter >= fi + lenr || end_iter < start_iter
              #println("oops")
              break
            end
            start_ts = time_ns()
            #println("Thread id: ", Threads.threadid(), " [", start_iter, ", ", end_iter, "], ", fi + lenr)
            for i in start_iter:end_iter
              local $(esc(index_variable)) = @inbounds i
              $(esc(body))
            end
            end_ts = time_ns()
            push!(tasks_log, LogInfo(start_ts, end_ts, Threads.threadid(), start_iter, end_iter))
          end
          outfile_name = "/Users/guillo0001/ghq/github.com/GuilloteauQ/LB.jl/lb4jl_thread_$(Threads.threadid()).csv" 
          open(outfile_name, "a") do file
            for task in tasks_log
              write(file, get_str(task, sched_v))
            end
          end
        end
      end
      end
      threading_run(thread_f, true)
    end
end

end
