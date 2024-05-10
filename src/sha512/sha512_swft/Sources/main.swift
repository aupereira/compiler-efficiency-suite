import Foundation

let BUFFER_SIZE: Int = 1024

let K: [UInt64] = [
    0x428a2f98d728ae22, 0x7137449123ef65cd, 0xb5c0fbcfec4d3b2f, 0xe9b5dba58189dbbc,
    0x3956c25bf348b538, 0x59f111f1b605d019, 0x923f82a4af194f9b, 0xab1c5ed5da6d8118,
    0xd807aa98a3030242, 0x12835b0145706fbe, 0x243185be4ee4b28c, 0x550c7dc3d5ffb4e2,
    0x72be5d74f27b896f, 0x80deb1fe3b1696b1, 0x9bdc06a725c71235, 0xc19bf174cf692694,
    0xe49b69c19ef14ad2, 0xefbe4786384f25e3, 0x0fc19dc68b8cd5b5, 0x240ca1cc77ac9c65,
    0x2de92c6f592b0275, 0x4a7484aa6ea6e483, 0x5cb0a9dcbd41fbd4, 0x76f988da831153b5,
    0x983e5152ee66dfab, 0xa831c66d2db43210, 0xb00327c898fb213f, 0xbf597fc7beef0ee4,
    0xc6e00bf33da88fc2, 0xd5a79147930aa725, 0x06ca6351e003826f, 0x142929670a0e6e70,
    0x27b70a8546d22ffc, 0x2e1b21385c26c926, 0x4d2c6dfc5ac42aed, 0x53380d139d95b3df,
    0x650a73548baf63de, 0x766a0abb3c77b2a8, 0x81c2c92e47edaee6, 0x92722c851482353b,
    0xa2bfe8a14cf10364, 0xa81a664bbc423001, 0xc24b8b70d0f89791, 0xc76c51a30654be30,
    0xd192e819d6ef5218, 0xd69906245565a910, 0xf40e35855771202a, 0x106aa07032bbd1b8,
    0x19a4c116b8d2d0c8, 0x1e376c085141ab53, 0x2748774cdf8eeb99, 0x34b0bcb5e19b48a8,
    0x391c0cb3c5c95a63, 0x4ed8aa4ae3418acb, 0x5b9cca4f7763e373, 0x682e6ff3d6b2b8a3,
    0x748f82ee5defb2fc, 0x78a5636f43172f60, 0x84c87814a1f0ab72, 0x8cc702081a6439ec,
    0x90befffa23631e28, 0xa4506cebde82bde9, 0xbef9a3f7b2c67915, 0xc67178f2e372532b,
    0xca273eceea26619c, 0xd186b8c721c0c207, 0xeada7dd6cde0eb1e, 0xf57d4f7fee6ed178,
    0x06f067aa72176fba, 0x0a637dc5a2c898a6, 0x113f9804bef90dae, 0x1b710b35131c471b,
    0x28db77f523047d84, 0x32caab7b40c72493, 0x3c9ebe0a15c9bebc, 0x431d67c49c100d4c,
    0x4cc5d4becb3e42b6, 0x597f299cfc657e2a, 0x5fcb6fab3ad6faec, 0x6c44198c4a475817,
]

func bytesToBigEndian(_ bytes: ArraySlice<UInt8>, _ index: Int) -> UInt64 {
    let part1 = UInt64(bytes[index]) << 56
    let part2 = UInt64(bytes[index + 1]) << 48
    let part3 = UInt64(bytes[index + 2]) << 40
    let part4 = UInt64(bytes[index + 3]) << 32
    let part5 = UInt64(bytes[index + 4]) << 24
    let part6 = UInt64(bytes[index + 5]) << 16
    let part7 = UInt64(bytes[index + 6]) << 8
    let part8 = UInt64(bytes[index + 7])
    
    return part1 | part2 | part3 | part4 | part5 | part6 | part7 | part8
}

func rotr(_ x: UInt64, _ n: UInt8) -> UInt64 {
    return (x >> n) | (x << (64 - n))
}

func ch(_ x: UInt64, _ y: UInt64, _ z: UInt64) -> UInt64 {
    return (x & y) ^ (~x & z)
}

func maj(_ x: UInt64, _ y: UInt64, _ z: UInt64) -> UInt64 {
    return (x & y) ^ (x & z) ^ (y & z)
}

func Sig0(_ x: UInt64) -> UInt64 {
    return rotr(x, 28) ^ rotr(x, 34) ^ rotr(x, 39)
}

func Sig1(_ x: UInt64) -> UInt64 {
    return rotr(x, 14) ^ rotr(x, 18) ^ rotr(x, 41)
}

func sig0(_ x: UInt64) -> UInt64 {
    return rotr(x, 1) ^ rotr(x, 8) ^ (x >> 7)
}

func sig1(_ x: UInt64) -> UInt64 {
    return rotr(x, (UInt8)(19)) ^ rotr(x, (UInt8)(61)) ^ (x >> 6)
}

func finalBlock(_ buffer: inout [UInt64], _ bytes: Int, _ totalBits: UInt64) {
    let block: Int = (bytes - 1) / 128
    var blockIndex: Int = block * 16
    let cByte: Int = bytes % 128
    let cWord: Int = cByte / 8

    if (cByte < 112) && (cByte > 0) {
        buffer[blockIndex + cWord] = buffer[blockIndex + cWord] & ~( 0xFF >> (8 * (cByte % 8))) | (0x80 << (8 * (7 - (cByte % 8))))
        for i in (cWord + 1)..<14 {
            buffer[blockIndex + i] = 0
        }
        buffer[blockIndex + 14] = 0
        buffer[blockIndex + 15] = totalBits
    }

    else if (cByte < 128) && (cByte > 0) {
        buffer[blockIndex + cWord] = buffer[blockIndex + cWord] & ~( 0xFF >> (8 * (cByte % 8))) | (0x80 << (8 * (7 - (cByte % 8))))
        if (cWord + 1) == 15 {
            buffer[15] = 0
        }

        if (block + 1 < BUFFER_SIZE) {
            blockIndex += 16
            for i in 0..<14 {
                buffer[blockIndex + i] = 0
            }
            buffer[blockIndex + 14] = 0
            buffer[blockIndex + 15] = totalBits
        }
    }

    else if (block + 1 < BUFFER_SIZE) {
        blockIndex += 16
        buffer[blockIndex] = 0x8000000000000000
        for i in 1..<14 {
            buffer[blockIndex + i] = 0
        }
        buffer[blockIndex + 14] = 0
        buffer[blockIndex + 15] = totalBits
    }
}

func addPadding(_ buffer: inout [UInt64], _ totalBits: UInt64, _ addEOS: Bool) {
    for i in 0..<14 {
        buffer[i] = 0
    }

    if addEOS {
        buffer[0] = 0x8000000000000000
    }

    buffer[14] = 0
    buffer[15] = totalBits
}

func sha512Rounds(_ buffer: inout [UInt64], _ state: inout [UInt64], _ iter: Int) {
    var a, b, c, d, e, f, g, h, t1, t2: UInt64

    var w: [UInt64] = [UInt64](repeating: 0, count: 80)

    var baseIndex: Int

    for block in 0..<iter {
        baseIndex = block * 16

        for i in 0..<16 {
            w[i] = buffer[baseIndex + i]
        }

        for i in 16..<80 {
            w[i] = sig1(w[i - 2]) &+ w[i - 7] &+ sig0(w[i - 15]) &+ w[i - 16]
        }

        a = state[0]
        b = state[1]
        c = state[2]
        d = state[3]
        e = state[4]
        f = state[5]
        g = state[6]
        h = state[7]

        for i in 0..<80 {
            t1 = h &+ Sig1(e) &+ ch(e, f, g) &+ K[i] &+ w[i]
            t2 = Sig0(a) &+ maj(a, b, c)
            h = g
            g = f
            f = e
            e = d &+ t1
            d = c
            c = b
            b = a
            a = t1 &+ t2
        }

        state[0] &+= a
        state[1] &+= b
        state[2] &+= c
        state[3] &+= d
        state[4] &+= e
        state[5] &+= f
        state[6] &+= g
        state[7] &+= h
    }
}

func sha512(_ filePath: String, _ iter: Int) {
    let totalBytes: UInt64
    let totalBits: UInt64

    let f: FileHandle = FileHandle(forReadingAtPath: filePath)!
    
    do {
        totalBytes = try f.seekToEnd()
        totalBits = totalBytes * 8
        try f.seek(toOffset: 0)
    } catch {
        print("Error getting file size")
        exit(1)
    }

    var buffers: [UInt8]
    var wbuffers: [UInt64] = [UInt64](repeating: 0, count: BUFFER_SIZE * 16)

    var state: [UInt64] = [
        0x6a09e667f3bcc908, 0xbb67ae8584caa73b,
        0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1,
        0x510e527fade682d1, 0x9b05688c2b3e6c1f,
        0x1f83d9abfb41bd6b, 0x5be0cd19137e2179
    ]

    if totalBytes == 0 {
        wbuffers[0] = 0x8000000000000000
        sha512Rounds(&wbuffers, &state, 1)
        for i in 0..<8 {
            print(String(format: "%016llx", state[i]), terminator: "")
        }
        exit(0);
    }

    for _ in 0..<((Int)((totalBytes - 1) / (UInt64)(BUFFER_SIZE * 128))) {
        buffers = [UInt8](f.readData(ofLength: BUFFER_SIZE * 128))
        for _ in 0..<iter {
            for j in 0..<BUFFER_SIZE * 16 {
                wbuffers[j] = bytesToBigEndian(buffers[j * 8..<j * 8 + 8], j * 8)
            }
            sha512Rounds(&wbuffers, &state, BUFFER_SIZE)
        }
    }

    var bytes: Int = (Int)((totalBytes) % (UInt64)(BUFFER_SIZE * 128))
    if bytes == 0 {
        bytes = BUFFER_SIZE * 128
    }
    let paddedWords = ((bytes - 1) / 8) + 1
    buffers = [UInt8](f.readData(ofLength: BUFFER_SIZE * 128))

    var paddedBuffers: [UInt8] = [UInt8](repeating: 0, count: paddedWords * 8)
    paddedBuffers[0..<bytes] = buffers[0..<bytes]
    for i in 0..<paddedWords {
        wbuffers[i] = bytesToBigEndian(paddedBuffers[i * 8..<i * 8 + 8], i * 8)
    }

    finalBlock(&wbuffers, bytes, totalBits)

    if (bytes > BUFFER_SIZE * 128 - 17) || (bytes == 0) {
        sha512Rounds(&wbuffers, &state, BUFFER_SIZE)
        addPadding(&wbuffers, totalBits, bytes % 128 == 0)
        sha512Rounds(&wbuffers, &state, 1)
    }
    else {
        sha512Rounds(&wbuffers, &state, (bytes + 144) / 128)
    }

    for i in 0..<8 {
        print(String(format: "%016llx", state[i]), terminator: "")
    }
    
    do {
        try f.close()
    } catch {
        print("Error closing file")
        exit(1)
    }
}

func main() {
    if CommandLine.arguments.count != 3{
        print("Usage: ./sha512_swft <File Path> <Iterations>")
        exit(1)
    }

    sha512(CommandLine.arguments[1], Int(CommandLine.arguments[2])!)
}

main()
