import Foundation

func fft(_ x: inout [Double]) {
    var real, imag, tReal, tImag: Double
    
    let n = x.count

    if n == 2 {
        return
    }

    var even: [Double] = Array(repeating: 0.0, count: n / 2)
    var odd: [Double] = Array(repeating: 0.0, count: n / 2)

    for i in stride(from: 0, to: n / 2, by: 2) {
        even[i] = x[i * 2]
        even[i + 1] = x[i * 2 + 1]
        odd[i] = x[i * 2 + 2]
        odd[i + 1] = x[i * 2 + 3]
    }

    fft(&even)
    fft(&odd)

    for k in stride(from: 0, to: n / 2, by: 2) {
        imag = -2 * Double.pi * Double(k) / Double(n)
        real = cos(imag)
        imag = sin(imag)

        tReal = real * odd[k] - imag * odd[k + 1]
        tImag = real * odd[k + 1] + imag * odd[k]

        x[k] = even[k] + tReal
        x[k + 1] = even[k + 1] + tImag

        x[k + n / 2] = even[k] - tReal
        x[k + n / 2 + 1] = even[k + 1] - tImag
    }
}

func fftLoop(size: Int, loops: Int) async {
    var x: [Double] = Array(repeating: 0.0, count: 2 * size)

    for _ in 0 ..< loops {
        for i in 0 ..< x.count {
            x[i] = Double.random(in: 0 ... 1)
        }
        fft(&x)
    }
}

func main() async {
    if CommandLine.arguments.count != 4 {
        print("Usage: ./fft_swft <Input Size> <Loops> <Threads>")
        return
    }

    let fftSize = 1 << Int(CommandLine.arguments[1])!
    let numLoops = Int(CommandLine.arguments[2])!
    let numThreads = Int(CommandLine.arguments[3])!

    if numThreads == 1 {
        await fftLoop(size: fftSize, loops: numLoops)
        exit(0)
    }

    var tasks = [Task<Void, Never>]()
    for _ in 0 ..< numThreads {
        let task = Task {
            await fftLoop(size: fftSize, loops: numLoops)
        }
        tasks.append(task)
    }

    for task in tasks {
        await task.value
    }
}

let semaphore = DispatchSemaphore(value: 0)

Task.detached {
    await main()
    semaphore.signal()
}

semaphore.wait()
