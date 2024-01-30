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
  tss
  fac2
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

function get_chunk(first_index::Int64, max_len::Int64, sched::Schemes)
  p = Threads.threadpoolsize()
  if sched == static::Schemes
    ((i::Int64) -> get_chunk_static(i, first_index, max_len, div(max_len, p)), p)
  elseif sched == dynamic::Schemes
    ((i::Int64) -> get_chunk_dynamic(i, first_index, max_len), max_len)
  elseif sched == gss::Schemes
    ((i::Int64) -> get_chunk_gss(i, first_index, max_len, p), ceil(log(p / max_len) / log((p-1)/p) + 1))
  elseif sched == tss::Schemes
    f = ceil(Int64, max_len / (2 * Threads.threadpoolsize()))
    l = 1
    S = ceil(Int64, 2 * max_len / (f + l))
    delta = floor(Int64, (f - l) / (S - 1))
    ((i::Int64) -> get_chunk_tss(i, first_index, max_len, f, l, S, delta), S)
  elseif sched == fac2::Schemes
    ((i::Int64) -> get_chunk_fac2(i, first_index, max_len, p), max_len)
  else
    println("Unknown sched: ", sched)
    throw("Unknown `sched`")
  end
end

macro mythreads(args...)
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
        # chunk_f, max_chunks = get_chunk(fi, lenr, gss::Schemes)
        chunk_f, max_chunks = get_chunk(fi, lenr, sched_v)
        queue_index = Threads.Atomic{Int}(1)

        function thread_f()
          while true
            index = Threads.atomic_add!(queue_index, 1)
            if index > max_chunks
              return
            end
            start_iter, end_iter = chunk_f(index)
            #end_iter = min(end_iter, fi + lenr)
            if start_iter >= fi + lenr || end_iter >= fi + lenr || end_iter < start_iter
              #println("oops")
              return
            end
            #println("Thread id: ", Threads.threadid(), " [", start_iter, ", ", end_iter, "], ", fi + lenr)
            for i in start_iter:end_iter
              local $(esc(index_variable)) = @inbounds i
              $(esc(body))
            end
          end
        end
      end
      end
      threading_run(thread_f, true)
    end
end


#function par_compute_pi(nb_iteration)::Float64
#  counts = Array{Bool}(undef, nb_iteration)
#  # @Threads.threads :static for i in 1:nb_iteration
#  @mythreads :dynamic for i in 1:nb_iteration
#    (x, y) = random_point()
#    counts[i] = (x * x + y * y < 1)
#  end
#  4 * sum(counts) / nb_iteration
#end



n = 41 * 1000000
#for i in 1:10 
#  pi_approx = @time compute_pi(n)
#  #println(pi_approx)
#end
#
#for i in 1:10
#  pi_approx_par =  @time par_compute_pi(n)
#  #println(pi_approx_par)
#end

function diverges(c::Complex, nb_iterations::Int64, epsilon::Float64)::Bool
  z = 0 + im * 0
  i = 0
  while abs(z) < epsilon && i < nb_iterations
    z = z * z + c
    i = i + 1
  end
  i < nb_iterations
end

function mandelbrot()
  xmin = -4.0
  xmax =  4.0
  ymin = -1.0
  ymax =  1.0
  delta = 0.001


  nb_rows = floor(Int64, (ymax - ymin) / delta)
  nb_cols = floor(Int64, (xmax - xmin) / delta)

  data = Matrix{Bool}(undef, nb_rows, nb_cols)

  for x in 1:nb_cols
    for y in 1:nb_rows
      c = ((xmax - xmin) * x / nb_cols + xmin) + im * ((ymax - ymin) * y / nb_rows + ymin)
      data[y, x] = diverges(c, 10000, 100.0)
    end
  end
  data
end

function mandelbrot_par(scheme::Schemes)
  xmin = -4.0
  xmax =  4.0
  ymin = -1.0
  ymax =  1.0
  delta = 0.001


  nb_rows = floor(Int64, (ymax - ymin) / delta)
  nb_cols = floor(Int64, (xmax - xmin) / delta)

  data = Matrix{Bool}(undef, nb_rows, nb_cols)

  @mythreads scheme for x in 1:nb_cols
    for y in 1:nb_rows
      c = ((xmax - xmin) * x / nb_cols + xmin) + im * ((ymax - ymin) * y / nb_rows + ymin)
      data[y, x] = diverges(c, 10000, 100.0)
    end
  end
  data
end

println("Mandelbrot --------------------")
print("Static ")
@time mandelbrot_par(static::Schemes)
print("GSS ")
@time mandelbrot_par(gss::Schemes)
print("TSS ")
@time mandelbrot_par(tss::Schemes)
print("Dynamic ")
@time mandelbrot_par(dynamic::Schemes)
print("FAC2 ")
@time mandelbrot_par(fac2::Schemes)



println("PI --------------------")

function par_compute_pi(nb_iteration, scheme::Schemes)::Float64
  counts = Array{Bool}(undef, nb_iteration)
  @mythreads scheme for i in 1:nb_iteration
    (x, y) = random_point()
    counts[i] = (x * x + y * y < 1)
  end
  4 * sum(counts) / nb_iteration
end

n = 10000000
print("Static ")
@time par_compute_pi(n, static::Schemes)
print("GSS ")
@time par_compute_pi(n, gss::Schemes)
print("TSS ")
@time par_compute_pi(n, tss::Schemes)
print("Dynamic ")
@time par_compute_pi(n, dynamic::Schemes)
print("FAC2 ")
@time par_compute_pi(n, fac2::Schemes)
