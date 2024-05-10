def rand_symmetric_ip(n, a)
    for i in 0...n
        a[i][i] = rand
        for j in (i+1)...n
            a[i][j] = rand
            a[j][i] = a[i][j]
        end
    end
end

def sgn(x)
    if x >= 0.0
        return 1.0
    else
        return -1.0
    end
end

def norm2(len, v)
    norm = 0.0

    for i in 0...len
        norm += v[i] * v[i]
    end

    return Math.sqrt(norm)
end

def v_div_s(len, v, s)
    for i in 0...len
        v[i] /= s
    end

    return v
end

def s2opai(len, v, a)
    for i in 0...len
        for j in 0...len
            a[i][j] = -2.0 * v[i] * v[j]
        end
    end

    for i in 0...len
        a[i][i] += 1.0
    end
end

def matrix_multiply(len, a, b, c, btmp)
    for i in 0...len
        for j in 0...len
            btmp[i][j] = b[j][i]
        end
    end

    for i in 0...len
        for j in 0...len
            sum = 0.0
            for k in 0...len
                sum += a[i][k] * btmp[j][k]
            end
            c[i][j] = sum
        end
    end
end

def tridiagonalize(n, a, p, temp, tpose, v)
    for k in 0...n-2
        for j in 0...k+1
            v[j] = 0.0
        end
        for j in k+1...n
            v[j] = a[j][k]
        end

        alpha = -sgn(v[k + 1]) * norm2(n, v)
        r = 2.0 * Math.sqrt(0.5 * ((alpha * alpha) - (alpha * v[k + 1])))

        v[k + 1] = (v[k + 1] - alpha) / r
        for j in k+2...n
            v[j] /= r
        end

        s2opai(n, v, p)
        matrix_multiply(n, p, a, temp, tpose)
        matrix_multiply(n, temp, p, a, tpose)
    end
end

def householder(len, a, a1, q, h, temp, tpose, v)
    a1 = a.map(&:dup)

    s = a1[len - 1][len - 1]
    t = a1[len - 2][len - 2]
    x = a1[len - 2][len - 1]
    d = (t - s) / 2.0
    shift = s - (x * x) / (d + sgn(d) * Math.sqrt(d * d + x * x))

    for i in 0...len
        a1[i][i] -= shift
    end

    for k in 0...len-1
        for j in 0...k
            v[j] = 0.0
        end
        for j in k...len
            v[j] = a1[j][k]
        end

        v[k] += norm2(len, v) * sgn(v[k + 1])
        v_div_s(len, v, norm2(len, v))

        if k == 0
            s2opai(len, v, q)
            matrix_multiply(len, q, a1, temp, tpose)
            a1, temp = temp, a1
        else
            s2opai(len, v, h)
            matrix_multiply(len, h, a1, temp, tpose)
            a1, temp = temp, a1
            matrix_multiply(len, q, h, temp, tpose)
            q, temp = temp, q
        end
    end

    matrix_multiply(len, a1, q, a, tpose)

    for i in 0...len
        a[i][i] += shift
    end
end

def eig_qr(iter, a, eig_vals, t1, t2, t3, t4, t5, v)
    ind = a.length
    
    while ind > 1
        for _ in 0...iter
            householder(ind, a, t1, t2, t3, t4, t5, v)
        end
        eig_vals[ind - 1] = a[ind - 1][ind - 1]
        ind -= 1
    end

    eig_vals[0] = a[0][0]

end

def eigen_loop(size, iter, loops)
    a = Array.new(size) { Array.new(size) }

    t1 = Array.new(size) { Array.new(size) }
    t2 = Array.new(size) { Array.new(size) }
    t3 = Array.new(size) { Array.new(size) }
    t4 = Array.new(size) { Array.new(size) }
    t5 = Array.new(size) { Array.new(size) }
    v = Array.new(size)

    eig_vals = Array.new(size)

    loops.times do
        rand_symmetric_ip(size, a)
        tridiagonalize(size, a, t1, t2, t3, v)
        eig_qr(iter, a, eig_vals, t1, t2, t3, t4, t5, v)
    end
end

def main
    if ARGV.length != 4
        puts "Usage: ruby eigenvalue_rb.rb <Matrix Size> <Convergence Iterations> <Loops> <Threads>"
        exit 1
    end

    m_size = ARGV[0].to_i
    c_iter = ARGV[1].to_i
    num_loops = ARGV[2].to_i
    num_threads = ARGV[3].to_i

    if num_threads == 1
        eigen_loop(m_size, c_iter, num_loops)
        exit 0
    end

    num_threads.times do |i|
        Process.spawn(RUBY_ENGINE, '-e', "require './#{__FILE__}'; eigen_loop(#{m_size}, #{c_iter}, #{num_loops})")
	end

    Process.waitall
end

main if __FILE__ == $PROGRAM_NAME
