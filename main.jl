function random_point()
  (rand(), rand())
end

function compute_pi(nb_iteration)::Float64
  count = 0
  for i in 1:nb_iteration
    (x, y) = random_point()
    if x * x + y * y < 1
      count += 1
    end
  end
  4 * count / nb_iteration
end

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

@enum Schemes begin
  static = 1
  dynamic = 2
  gss
end

function get_chunk_dynamic(j::Int64, first_index::Int64, _max_len::Int64)::Tuple{Int, Int}
  (first_index + (j - 1), first_index + j)
end

function get_chunk_static(j::Int64, first_index::Int64, max_len::Int64)::Tuple{Int, Int}
  (floor(first_index + (j - 1) * max_len / Threads.threadpoolsize()), floor(first_index + j * max_len / Threads.threadpoolsize() - 1))
end

function get_chunk_gss(j::Int64, first_index::Int64, max_len::Int64)::Tuple{Int, Int}
  p = Threads.threadpoolsize()
  offset = max_len - max_len * ((p-1)/p)^(j-1)
  k = (max_len / p) * ((p-1)/p)^(j-1)
  (floor(first_index + offset), floor(first_index + offset + k - 1))
end

function get_chunk(first_index::Int64, max_len::Int64, sched::Schemes)
  p = Threads.threadpoolsize()
  if sched == static::Schemes
    ((i::Int64) -> get_chunk_static(i, first_index, max_len), p)
  elseif sched == dynamic::Schemes
    ((i::Int64) -> get_chunk_dynamic(i, first_index, max_len), max_len)
  elseif sched == gss::Schemes
    ((i::Int64) -> get_chunk_gss(i, first_index, max_len), ceil(log(p / max_len) / log((p-1)/p) + 1))
  else
    println("Unknown sched: ", sched)
    throw("Unknown `sched`")
  end
end

macro mythreads(args...)
    println(args)
    sched, loop = args
    println(sched)
    sched_v = sched.value
    iterator = loop.args[1]
    body = loop.args[2]
    index_variable = iterator.args[1]
    range = iterator.args[2]

    println(index_variable)

    quote
      local thread_f
      let range = $(esc(range))
        lenr = length(range)
        fi = firstindex(range)
        chunk_f, max_chunks = get_chunk(fi, lenr, static::Schemes)
        queue_index = Threads.Atomic{Int}(1)

        function thread_f()
          while true
            index = Threads.atomic_add!(queue_index, 1)
            if index > max_chunks
              return
            end
            start_iter, end_iter = chunk_f(index)
            end_iter = min(end_iter, fi + lenr)
            if end_iter < start_iter
              println("oops")
              return
            end
            #println("Thread id: ", Threads.threadid(), " [", start_iter, ", ", end_iter, "]")
            for i in start_iter:end_iter
              local $(esc(index_variable)) = @inbounds i
              $(esc(body))
            end
          end
        end
      end
      threading_run(thread_f, true)
    end
end


function par_compute_pi(nb_iteration)::Float64
  counts = Array{Bool}(undef, nb_iteration)
  # @Threads.threads :static for i in 1:nb_iteration
  @mythreads :dynamic for i in 1:nb_iteration
    (x, y) = random_point()
    counts[i] = (x * x + y * y < 1)
  end
  4 * sum(counts) / nb_iteration
end



n = 3 * 13 * 1000000
for i in 1:10 
  pi_approx = @time compute_pi(n)
  println(pi_approx)
end

for i in 1:10
  pi_approx_par =  @time par_compute_pi(n)
  println(pi_approx_par)
end
