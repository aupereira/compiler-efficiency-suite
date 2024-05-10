import os
import sys

BUFFER_SIZE = 1024

K = [
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
    0x4cc5d4becb3e42b6, 0x597f299cfc657e2a, 0x5fcb6fab3ad6faec, 0x6c44198c4a475817
]


def to_big_endian(n):
    return (
        (n >> 56) | (n >> 40) & 0x000000000000FF00 |
		(n >> 24) & 0x0000000000FF0000 | (n >> 8) & 0x00000000FF000000 |
		(n << 8) & 0x000000FF00000000 | (n << 24) & 0x0000FF0000000000 |
		(n << 40) & 0x00FF000000000000 | (n << 56)
    ) & 0xFFFFFFFFFFFFFFFF


def rotr(x, n):
    return ((x >> n) | (x << (64 - n))) & 0xFFFFFFFFFFFFFFFF


def ch(x, y, z):
    return (x & y) ^ (~x & z)


def maj(x, y, z):
    return (x & y) ^ (x & z) ^ (y & z)


def Sig0(x):
    return rotr(x, 28) ^ rotr(x, 34) ^ rotr(x, 39)


def Sig1(x):
    return rotr(x, 14) ^ rotr(x, 18) ^ rotr(x, 41)


def sig0(x):
    return rotr(x, 1) ^ rotr(x, 8) ^ (x >> 7)


def sig1(x):
    return rotr(x, 19) ^ rotr(x, 61) ^ (x >> 6)


def final_block(buffer, bytes, total_bits):
    block = (bytes - 1) // 128
    block_i = block * 16
    current_byte = bytes % 128
    current_word = current_byte // 8

    if current_byte < 112 and current_byte > 0:
        buffer[block_i + current_word] = (buffer[block_i + current_word] << 8 * (8 - (current_byte % 8))) & ~(0xFF >> (8 * (current_byte % 8))) | (0x80 << (8 * (7 - (current_byte % 8))))
        for i in range(current_word + 1, 14):
            buffer[block_i + i] = 0
        buffer[block_i + 14] = 0
        buffer[block_i + 15] = total_bits
    
    elif current_byte < 128 and current_byte > 0:
        buffer[block_i + current_word] = (buffer[block_i + current_word] << 8 * (8 - (current_byte % 8))) & ~(0xFF >> (8 * (current_byte % 8))) | (0x80 << (8 * (7 - (current_byte % 8))))
        if current_byte + 1 == 15:
            buffer[block_i + 15] = 0
        
        if block + 1 < BUFFER_SIZE:
            block_i += 16
            for i in range(block_i, block_i + 14):
                buffer[i] = 0
            buffer[block_i + 14] = 0
            buffer[block_i + 15] = total_bits
    
    elif block + 1 < BUFFER_SIZE:
        block_i += 16
        buffer[block_i] = 0x8000000000000000
        for i in range(block_i + 1, block_i + 14):
            buffer[i] = 0
        buffer[block_i + 14] = 0
        buffer[block_i + 15] = total_bits

    return buffer


def add_padding(buffer, total_bits, add_EOS):
    for i in range(0, 14):
        buffer[i] = 0

    if add_EOS:
        buffer[0] = 0x8000000000000000

    buffer[14] = 0
    buffer[15] = total_bits

    return buffer


def sha512_rounds(buffer, state, iter):
    a, b, c, d, e, f, g, h, t1, t2 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    w = [0] * 80
    
    for block in range(0, iter):
        base = block * 16

        for i in range(0, 16):
            w[i] = buffer[base + i]
        
        for i in range(16, 80):
            w[i] = (sig1(w[i - 2]) + w[i - 7] + sig0(w[i - 15]) + w[i - 16]) & 0xFFFFFFFFFFFFFFFF

        a = state[0]
        b = state[1]
        c = state[2]
        d = state[3]
        e = state[4]
        f = state[5]
        g = state[6]
        h = state[7]

        for i in range(0, 80):
            t1 = (h + Sig1(e) + ch(e, f, g) + K[i] + w[i]) & 0xFFFFFFFFFFFFFFFF
            t2 = (Sig0(a) + maj(a, b, c)) & 0xFFFFFFFFFFFFFFFF
            h = g
            g = f
            f = e
            e = (d + t1) & 0xFFFFFFFFFFFFFFFF
            d = c
            c = b
            b = a
            a = (t1 + t2) & 0xFFFFFFFFFFFFFFFF

        state[0] = (state[0] + a) & 0xFFFFFFFFFFFFFFFF
        state[1] = (state[1] + b) & 0xFFFFFFFFFFFFFFFF
        state[2] = (state[2] + c) & 0xFFFFFFFFFFFFFFFF
        state[3] = (state[3] + d) & 0xFFFFFFFFFFFFFFFF
        state[4] = (state[4] + e) & 0xFFFFFFFFFFFFFFFF
        state[5] = (state[5] + f) & 0xFFFFFFFFFFFFFFFF
        state[6] = (state[6] + g) & 0xFFFFFFFFFFFFFFFF
        state[7] = (state[7] + h) & 0xFFFFFFFFFFFFFFFF

    return state


def sha512(file_path, iter):
    qbuffers = [0] * (BUFFER_SIZE * 16)
    
    state = [
        0x6a09e667f3bcc908, 0xbb67ae8584caa73b,
        0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1,
        0x510e527fade682d1, 0x9b05688c2b3e6c1f,
        0x1f83d9abfb41bd6b, 0x5be0cd19137e2179
    ]

    try:
        f = open(file_path, "rb")
    except:
        print("Error opening file")
        os.exit(1)

    f.seek(0, os.SEEK_END)
    total_bytes = f.tell()
    total_bits = total_bytes * 8
    f.seek(0, os.SEEK_SET)

    if total_bytes == 0:
        qbuffers[0] = 0x8000000000000000
        sha512_rounds(qbuffers, state, 1)
        for i in range(0, 8):
            print(f'%016x' % state[i], end='')
        return

    for _ in range(0, (total_bytes - 1) // (BUFFER_SIZE * 128)):
        buffers = f.read(BUFFER_SIZE * 128)
        for j in range(iter):
            for i in range(BUFFER_SIZE * 16):
                qbuffers[i] = int.from_bytes(
                    buffers[8*i:8*i+8], byteorder='big', signed=False)
            state = sha512_rounds(qbuffers, state, BUFFER_SIZE)

    buffers = f.read(BUFFER_SIZE * 128)
    bytes = total_bytes % (BUFFER_SIZE * 128)
    for i in range(0, ((bytes - 1) // 8) + 1):
        qbuffers[i] = int.from_bytes(
            buffers[8*i:8*i+8], byteorder='big', signed=False)
    qbuffers = final_block(qbuffers, bytes, total_bits)
    
    if bytes > BUFFER_SIZE * 128 - 17:
        state = sha512_rounds(qbuffers, state, BUFFER_SIZE)
        qbuffers = add_padding(qbuffers, total_bits, bytes % 128 == 0)
        state = sha512_rounds(qbuffers, state, 1)
    else:
        state = sha512_rounds(qbuffers, state, (bytes + 144) // 128)

    for i in range(0, 8):
        print(f'%016x' % state[i], end='')

    try:
        f.close()
    except:
        print("Error closing file")
        os.exit(1)


def main():
    if len(sys.argv) != 3:
        print("Usage: py main.py <File Path> <Iterations>")
        os.exit(1)
    
    sha512(sys.argv[1], int(sys.argv[2]))


if __name__ == "__main__":
    main()
