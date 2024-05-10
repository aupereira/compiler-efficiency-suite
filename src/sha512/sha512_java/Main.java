import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;

class Main {
  private static final int BUFFER_SIZE = 1024;

  private static final long[] K = {
    0x428a2f98d728ae22L, 0x7137449123ef65cdL, 0xb5c0fbcfec4d3b2fL, 0xe9b5dba58189dbbcL,
    0x3956c25bf348b538L, 0x59f111f1b605d019L, 0x923f82a4af194f9bL, 0xab1c5ed5da6d8118L,
    0xd807aa98a3030242L, 0x12835b0145706fbeL, 0x243185be4ee4b28cL, 0x550c7dc3d5ffb4e2L,
    0x72be5d74f27b896fL, 0x80deb1fe3b1696b1L, 0x9bdc06a725c71235L, 0xc19bf174cf692694L,
    0xe49b69c19ef14ad2L, 0xefbe4786384f25e3L, 0x0fc19dc68b8cd5b5L, 0x240ca1cc77ac9c65L,
    0x2de92c6f592b0275L, 0x4a7484aa6ea6e483L, 0x5cb0a9dcbd41fbd4L, 0x76f988da831153b5L,
    0x983e5152ee66dfabL, 0xa831c66d2db43210L, 0xb00327c898fb213fL, 0xbf597fc7beef0ee4L,
    0xc6e00bf33da88fc2L, 0xd5a79147930aa725L, 0x06ca6351e003826fL, 0x142929670a0e6e70L,
    0x27b70a8546d22ffcL, 0x2e1b21385c26c926L, 0x4d2c6dfc5ac42aedL, 0x53380d139d95b3dfL,
    0x650a73548baf63deL, 0x766a0abb3c77b2a8L, 0x81c2c92e47edaee6L, 0x92722c851482353bL,
    0xa2bfe8a14cf10364L, 0xa81a664bbc423001L, 0xc24b8b70d0f89791L, 0xc76c51a30654be30L,
    0xd192e819d6ef5218L, 0xd69906245565a910L, 0xf40e35855771202aL, 0x106aa07032bbd1b8L,
    0x19a4c116b8d2d0c8L, 0x1e376c085141ab53L, 0x2748774cdf8eeb99L, 0x34b0bcb5e19b48a8L,
    0x391c0cb3c5c95a63L, 0x4ed8aa4ae3418acbL, 0x5b9cca4f7763e373L, 0x682e6ff3d6b2b8a3L,
    0x748f82ee5defb2fcL, 0x78a5636f43172f60L, 0x84c87814a1f0ab72L, 0x8cc702081a6439ecL,
    0x90befffa23631e28L, 0xa4506cebde82bde9L, 0xbef9a3f7b2c67915L, 0xc67178f2e372532bL,
    0xca273eceea26619cL, 0xd186b8c721c0c207L, 0xeada7dd6cde0eb1eL, 0xf57d4f7fee6ed178L,
    0x06f067aa72176fbaL, 0x0a637dc5a2c898a6L, 0x113f9804bef90daeL, 0x1b710b35131c471bL,
    0x28db77f523047d84L, 0x32caab7b40c72493L, 0x3c9ebe0a15c9bebcL, 0x431d67c49c100d4cL,
    0x4cc5d4becb3e42b6L, 0x597f299cfc657e2aL, 0x5fcb6fab3ad6faecL, 0x6c44198c4a475817L,
  };

  static long rotr(long x, byte n) {
    return (x >>> n) | (x << (64 - n));
  }

  static long ch(long x, long y, long z) {
    return (x & y) ^ (~x & z);
  }

  static long maj(long x, long y, long z) {
    return (x & y) ^ (x & z) ^ (y & z);
  }

  static long Sig0(long x) {
    return rotr(x, (byte) 28) ^ rotr(x, (byte) 34) ^ rotr(x, (byte) 39);
  }

  static long Sig1(long x) {
    return rotr(x, (byte) 14) ^ rotr(x, (byte) 18) ^ rotr(x, (byte) 41);
  }

  static long sig0(long x) {
    return rotr(x, (byte) 1) ^ rotr(x, (byte) 8) ^ (x >>> 7);
  }

  static long sig1(long x) {
    return rotr(x, (byte) 19) ^ rotr(x, (byte) 61) ^ (x >>> 6);
  }

  static void finalBlock(long[] buffer, int bytes, long totalBits) {
    int block = (bytes - 1) / 128;
    int blockIndex = block * 16;
    int cByte = bytes % 128;
    int cQword = cByte / 8;

    if ((cByte < 112) & (cByte > 0)) {
      buffer[blockIndex + cQword] =
          buffer[blockIndex + cQword] & ~(0xFFL >> (8 * (cByte % 8)))
              | 0x80L << (8 * (7 - (cByte % 8)));
      for (int i = cQword + 1; i < 16; i++) {
        buffer[blockIndex + i] = 0;
      }
      buffer[blockIndex + 14] = 0;
      buffer[blockIndex + 15] = totalBits;
    } else if ((cByte < 128) & (cByte > 0)) {
      buffer[blockIndex + cQword] =
          buffer[blockIndex + cQword] & ~(0xFFL >> (8 * (cByte % 8)))
              | 0x80L << (8 * (7 - (cByte % 8)));
      if (cQword + 1 == 15) {
        buffer[blockIndex + 15] = 0;
      }

      if (block < BUFFER_SIZE - 1) {
        blockIndex += 16;
        for (int i = blockIndex; i < 14; i++) {
          buffer[blockIndex + i] = 0;
        }
        buffer[blockIndex + 14] = 0;
        buffer[blockIndex + 15] = totalBits;
      }
    } else if (block < BUFFER_SIZE - 1) {
      blockIndex += 16;
      buffer[blockIndex] = 0x8000000000000000L;
      for (int i = 1; i < 14; i++) {
        buffer[blockIndex + i] = 0;
      }
      buffer[blockIndex + 14] = 0;
      buffer[blockIndex + 15] = totalBits;
    }
  }

  static void addPadding(long[] buffer, long totalBits, boolean addEOS) {
    for (int i = 0; i < 14; i++) {
      buffer[i] = 0;
    }

    if (addEOS) {
      buffer[0] = 0x8000000000000000L;
    }

    buffer[14] = 0;
    buffer[15] = totalBits;
  }

  static void sha512Rounds(long[] buffer, long[] state, int iter) {
    long a, b, c, d, e, f, g, h, t1, t2;
    long[] w = new long[80];
    int base;

    for (int block = 0; block < iter; block++) {
      base = block * 16;

      for (int i = 0; i < 16; i++) {
        w[i] = buffer[base + i];
      }

      for (int i = 16; i < 80; i++) {
        w[i] = sig1(w[i - 2]) + w[i - 7] + sig0(w[i - 15]) + w[i - 16];
      }

      a = state[0];
      b = state[1];
      c = state[2];
      d = state[3];
      e = state[4];
      f = state[5];
      g = state[6];
      h = state[7];

      for (int i = 0; i < 80; i++) {
        t1 = h + Sig1(e) + ch(e, f, g) + K[i] + w[i];
        t2 = Sig0(a) + maj(a, b, c);
        h = g;
        g = f;
        f = e;
        e = d + t1;
        d = c;
        c = b;
        b = a;
        a = t1 + t2;
      }

      state[0] += a;
      state[1] += b;
      state[2] += c;
      state[3] += d;
      state[4] += e;
      state[5] += f;
      state[6] += g;
      state[7] += h;
    }
  }

  static void sha512(String filePath, int iter) {
    int bytes = 0;
    long totalBytes = 0;
    long totalBits = 0;

    FileInputStream f = null;
    try {
      f = new FileInputStream(filePath);
    } catch (FileNotFoundException e) {
      System.out.println("Error opening file.");
      System.exit(1);
    } catch (SecurityException e) {
      System.out.println("Error opening file.");
      System.exit(1);
    }

    File fi = null;
    try {
      fi = new File(filePath);
      totalBytes = fi.length();
      totalBits = totalBytes * 8;
    } catch (NullPointerException e) {
      System.out.println("Error getting file size.");
      System.exit(1);
    } catch (SecurityException e) {
      System.out.println("Error getting file size.");
      System.exit(1);
    }

    byte[] buffers = new byte[BUFFER_SIZE * 128];
    long[] qbuffers = new long[BUFFER_SIZE * 16];
    ByteBuffer bb = ByteBuffer.wrap(buffers);

    long[] state = {
      0x6a09e667f3bcc908L, 0xbb67ae8584caa73bL,
      0x3c6ef372fe94f82bL, 0xa54ff53a5f1d36f1L,
      0x510e527fade682d1L, 0x9b05688c2b3e6c1fL,
      0x1f83d9abfb41bd6bL, 0x5be0cd19137e2179L
    };

    if (totalBytes == 0) {
      qbuffers[0] = 0x8000000000000000L;
      sha512Rounds(qbuffers, state, 1);
      for (int i = 0; i < 8; i++) {
        System.out.printf("%016x", state[i]);
      }
      System.exit(0);
    }

    for (int _i = 0; _i < (totalBytes - 1) / (BUFFER_SIZE * 128); _i++) {
      try {
        bytes = f.read(buffers);
      } catch (IOException e) {
        System.out.println("Error reading file.");
        System.exit(1);
      }
      for (int _j = 0; _j < iter; _j++) {
        bb.position(0);
        for (int i = 0; i < BUFFER_SIZE * 16; i++) {
          qbuffers[i] = bb.getLong();
        }
        sha512Rounds(qbuffers, state, BUFFER_SIZE);
      }
    }

    try {
      bytes = f.read(buffers);
    } catch (IOException e) {
      System.out.println("Error reading file.");
      System.exit(1);
    }
    bb.position(0);
    for (int i = 0; i < ((bytes - 1) / 8) + 1; i++) {
      qbuffers[i] = bb.getLong();
    }
    finalBlock(qbuffers, bytes, totalBits);

    if (bytes > BUFFER_SIZE * 128 - 17) {
      sha512Rounds(qbuffers, state, BUFFER_SIZE);
      addPadding(qbuffers, totalBits, bytes % 128 == 0);
      sha512Rounds(qbuffers, state, 1);
    } else {
      sha512Rounds(qbuffers, state, (bytes + 144) / 128);
    }

    for (int i = 0; i < 8; i++) {
      System.out.printf("%016x", state[i]);
    }

    try {
      f.close();
    } catch (IOException e) {
      System.out.println("Error closing file.");
      System.exit(1);
    }
  }

  public static void main(String[] args) {
    if (args.length != 2) {
      System.out.println("Usage: java Main <File Path> <Iterations>");
      System.exit(1);
    }

    sha512(args[0], Integer.parseInt(args[1]));
  }
}
