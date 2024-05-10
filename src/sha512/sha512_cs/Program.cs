using System.Buffers.Binary;
using System.IO;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

class Program
{

    const int BUFFER_SIZE = 1024;

    [System.Runtime.InteropServices.StructLayout(LayoutKind.Explicit)]
    unsafe struct Buffer
    {
        [System.Runtime.InteropServices.FieldOffset(0)]
        public fixed byte bytes[BUFFER_SIZE * 128];
        [System.Runtime.InteropServices.FieldOffset(0)]
        public fixed ulong words[BUFFER_SIZE * 16];
    }

    static ulong ROTR(ulong x, byte n)
    {
        return (x >> n) | (x << (64 - n));
    }

    static ulong CH(ulong x, ulong y, ulong z)
    {
        return (x & y) ^ (~x & z);
    }

    static ulong MAJ(ulong x, ulong y, ulong z)
    {
        return (x & y) ^ (x & z) ^ (y & z);
    }

    static ulong SIG0(ulong x)
    {
        return ROTR(x, 28) ^ ROTR(x, 34) ^ ROTR(x, 39);
    }

    static ulong SIG1(ulong x)
    {
        return ROTR(x, 14) ^ ROTR(x, 18) ^ ROTR(x, 41);
    }

    static ulong sig0(ulong x)
    {
        return ROTR(x, 1) ^ ROTR(x, 8) ^ (x >> 7);
    }

    static ulong sig1(ulong x)
    {
        return ROTR(x, 19) ^ ROTR(x, 61) ^ (x >> 6);
    }

    unsafe static void FinalBlock(ref Buffer buffer, int bytes, ulong totalBits)
    {
        int block = (bytes - 1) / 128;
        int blockIndex = block * 16;
        int cByte = bytes % 128;
        int cWord = cByte / 8;

        if ((cByte < 112) & (cByte > 0))
        {
            buffer.bytes[cByte] = 0x80;
            for (int i = cByte + 1; i < 112; i++)
            {
                buffer.bytes[i] = 0;
            }
            buffer.words[blockIndex + 14] = 0;
            buffer.words[blockIndex + 15] = BinaryPrimitives.ReverseEndianness(totalBits);
        }

        else if ((cByte < 128) & (cByte > 0))
        {
            buffer.bytes[cByte] = 0x80;
            for (int i = cByte + 1; i < 128; i++)
            {
                buffer.bytes[i] = 0;
            }

            if (block < BUFFER_SIZE - 1)
            {
                blockIndex += 16;
                for (int i = 0; i < 14; i++)
                {
                    buffer.words[blockIndex + i] = 0;
                }
                buffer.words[blockIndex + 14] = 0;
                buffer.words[blockIndex + 15] = BinaryPrimitives.ReverseEndianness(totalBits);
            }
        }

        else if (block < BUFFER_SIZE - 1)
        {
            blockIndex += 16;
            for (int i = 0; i < 14; i++)
            {
                buffer.words[blockIndex + i] = 0;
            }
            buffer.bytes[blockIndex * 8] = 0x80;
            buffer.words[blockIndex + 14] = 0;
            buffer.words[blockIndex + 15] = BinaryPrimitives.ReverseEndianness(totalBits);
        }
    }

    unsafe static void AddPadding(ref Buffer buffer, ulong totalBits, bool addEOS)
    {
        for (int i = 0; i < 14; i++)
        {
            buffer.words[i] = 0;
        }

        if (addEOS)
        {
            buffer.bytes[0] = 0x80;
        }

        buffer.words[14] = 0;
        buffer.words[15] = BinaryPrimitives.ReverseEndianness(totalBits);
    }

    unsafe static void SHA512Rounds(ref Buffer buffer, ulong[] state, int iter)
    {
        ulong[] K = {
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
        };

        ulong a, b, c, d, e, f, g, h, t1, t2;
        ulong[] w = new ulong[80];

        int baseIndex;

        for (int block = 0; block < iter; block++)
        {
            baseIndex = block * 16;

            for (int i = 0; i < 16; i++)
            {
                w[i] = BinaryPrimitives.ReverseEndianness(buffer.words[baseIndex + i]);
            }

            for (int i = 16; i < 80; i++)
            {
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

            for (int i = 0; i < 80; i++)
            {
                t1 = h + SIG1(e) + CH(e, f, g) + K[i] + w[i];
                t2 = SIG0(a) + MAJ(a, b, c);
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

    unsafe static void SHA512(string filePath, int iter)
    {
        FileStream fs = File.OpenRead(filePath);

        long totalBytes = fs.Length;
        ulong totalBits = (ulong)totalBytes * 8;

        Buffer buffers = new Buffer();
        Span<byte> span = new Span<byte>(Unsafe.AsPointer(ref buffers.bytes[0]), BUFFER_SIZE * 128);

        ulong[] state = new ulong[8] {
            0x6a09e667f3bcc908, 0xbb67ae8584caa73b,
            0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1,
            0x510e527fade682d1, 0x9b05688c2b3e6c1f,
            0x1f83d9abfb41bd6b, 0x5be0cd19137e2179
        };

        if (totalBytes == 0)
        {
            buffers.bytes[0] = 0x80;
            SHA512Rounds(ref buffers, state, 1);
            for (int i = 0; i < 8; i++)
            {
                Console.Write("{0:x16}", state[i]);
            }
            Environment.Exit(0);
        }

        for (int i = 0; i < (totalBytes - 1) / (BUFFER_SIZE * 128); i++)
        {
            fs.ReadExactly(span);
            for (int j = 0; j < iter; j++)
            {
                SHA512Rounds(ref buffers, state, BUFFER_SIZE);
            }
        }

        int bytes = (int)(totalBytes % (BUFFER_SIZE * 128));
        span = new Span<byte>(Unsafe.AsPointer(ref buffers.bytes[0]), bytes);
        fs.ReadExactly(span);
        FinalBlock(ref buffers, bytes, totalBits);

        if (bytes > BUFFER_SIZE * 128 - 17)
        {
            SHA512Rounds(ref buffers, state, BUFFER_SIZE);
            AddPadding(ref buffers, totalBits, bytes % 128 == 0);
            SHA512Rounds(ref buffers, state, 1);
        }
        else
        {
            SHA512Rounds(ref buffers, state, (bytes + 144) / 128);
        }

        for (int i = 0; i < 8; i++)
        {
            Console.Write("{0:x16}", state[i]);
        }

        fs.Close();
    }

    static void Main(string[] args)
    {
        if (args.Length != 2)
        {
            Console.WriteLine("Usage: ./sha512_cs <File Path> <Iterations>");
            return;
        }

        SHA512(args[0], int.Parse(args[1]));
    }
}
