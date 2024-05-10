# frozen_string_literal: true

BUFFER_SIZE = 6144

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
].freeze

def rotr(x, n)
  ((x >> n) | (x << (64 - n)))
end

def ch(x, y, z)
  (x & y) ^ (~x & z)
end

def maj(x, y, z)
  (x & y) ^ (x & z) ^ (y & z)
end

def Sig0(x)
  rotr(x, 28) ^ rotr(x, 34) ^ rotr(x, 39)
end

def Sig1(x)
  rotr(x, 14) ^ rotr(x, 18) ^ rotr(x, 41)
end

def sig0(x)
  rotr(x, 1) ^ rotr(x, 8) ^ (x >> 7)
end

def sig1(x)
  rotr(x, 19) ^ rotr(x, 61) ^ (x >> 6)
end

def final_block(buffer, bytes, total_bits)
  block = (bytes - 1) / 128
  block_i = block * 16
  c_byte = bytes % 128
  c_word = c_byte / 8

  if c_byte < 112 && c_byte != 0
    buffer[block_i + c_word] ||= 0
    buffer[block_i + c_word] = buffer[block_i + c_word] | (0x8000000000000000 >> (8 * (c_byte % 8)))
    (c_word + 1...14).each do |i|
      buffer[block_i + i] = 0
    end
    buffer[block_i + 14] = 0
    buffer[block_i + 15] = total_bits

  elsif c_byte != 0
    buffer[block_i + c_word] ||= 0
    buffer[block_i + c_word] = buffer[block_i + c_word] | (0x8000000000000000 >> (8 * (c_byte % 8)))
    buffer[block_i + 15] = 0 if c_word == 14

    if block < BUFFER_SIZE - 1
      block_i += 16
      (block_i...block_i + 14).each do |i|
        buffer[i] = 0
      end
      buffer[block_i + 14] = 0
      buffer[block_i + 15] = total_bits
    end

  elsif block < BUFFER_SIZE - 1
    block_i += 16
    buffer[block_i] = 0x8000000000000000
    (block_i + 1...block_i + 14).each do |i|
      buffer[i] = 0
    end
    buffer[block_i + 14] = 0
    buffer[block_i + 15] = total_bits
  end
end

def add_padding(buffer, total_bits, add_eos)
  (0...14).each do |i|
    buffer[i] = 0
  end

  buffer[0] = 0x8000000000000000 if add_eos

  buffer[14] = 0
  buffer[15] = total_bits
end

def sha512_rounds(buffer, state, iter)
  w = Array.new(80)

  (0...iter).each do |block|
    base = block * 16

    (0...16).each do |i|
      w[i] = buffer[base + i]
    end

    (16...80).each do |i|
      w[i] = (sig1(w[i - 2]) + w[i - 7] + sig0(w[i - 15]) + w[i - 16]) & 0xFFFFFFFFFFFFFFFF
    end

    a = state[0]
    b = state[1]
    c = state[2]
    d = state[3]
    e = state[4]
    f = state[5]
    g = state[6]
    h = state[7]

    (0...80).each do |i|
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
    end

    state[0] = (state[0] + a) & 0xFFFFFFFFFFFFFFFF
    state[1] = (state[1] + b) & 0xFFFFFFFFFFFFFFFF
    state[2] = (state[2] + c) & 0xFFFFFFFFFFFFFFFF
    state[3] = (state[3] + d) & 0xFFFFFFFFFFFFFFFF
    state[4] = (state[4] + e) & 0xFFFFFFFFFFFFFFFF
    state[5] = (state[5] + f) & 0xFFFFFFFFFFFFFFFF
    state[6] = (state[6] + g) & 0xFFFFFFFFFFFFFFFF
    state[7] = (state[7] + h) & 0xFFFFFFFFFFFFFFFF
  end
end

def sha512(file_path, iter)
  f = File.open(file_path, 'rb')
  total_bytes = File.size(file_path)

  qbuffers = Array.new(BUFFER_SIZE * 16, 0)

  state = [
    0x6a09e667f3bcc908, 0xbb67ae8584caa73b,
    0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1,
    0x510e527fade682d1, 0x9b05688c2b3e6c1f,
    0x1f83d9abfb41bd6b, 0x5be0cd19137e2179
  ]

  if total_bytes.zero?
    qbuffers[0] = 0x8000000000000000
    sha512_rounds(qbuffers, state, 1)
    state.each { |x| print format('%016x', x) }
    exit 0
  end

  total_bits = total_bytes * 8

  (0...(total_bytes - 1) / (BUFFER_SIZE * 128)).each do
    buffers = f.read(BUFFER_SIZE * 128)
    (0...iter).each do
      qbuffers.replace(buffers.unpack('Q>*'))
      sha512_rounds(qbuffers, state, BUFFER_SIZE)
    end
  end

  bytes = total_bytes % (BUFFER_SIZE * 128)
  bytes = BUFFER_SIZE * 128 if bytes.zero?

  buffers = f.read(bytes)
  if bytes % 8 != 0
    (bytes...bytes + (8 - (bytes % 8))).each do |i|
      buffers[i] = "\0"
    end
  end
  (0...((bytes - 1) / 8) + 1).each do |i|
    qbuffers[i] = buffers.slice(i * 8, i * 8 + 8).unpack1('Q>')
  end
  final_block(qbuffers, bytes, total_bits)

  if bytes > BUFFER_SIZE * 128 - 17
    sha512_rounds(qbuffers, state, BUFFER_SIZE)
    add_padding(qbuffers, total_bits, (bytes % 128).zero?)
    sha512_rounds(qbuffers, state, 1)
  else
    sha512_rounds(qbuffers, state, (bytes + 144) / 128)
  end

  state.each { |x| print format('%016x', x) }

  f.close
end

def main
  if ARGV.length != 2
    puts 'Usage: ruby main.rb <File Path> <Iterations>'
    exit 1
  end

  sha512(ARGV[0], ARGV[1].to_i)
end

main if __FILE__ == $PROGRAM_NAME
