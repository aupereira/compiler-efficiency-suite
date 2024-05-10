use std::env;
use std::fs::File;
use std::io::prelude::*;

const BUFFER_SIZE: usize = 1024;

#[repr(C)]
#[derive(Copy, Clone)]
union Block {
    bytes: [u8; 128],
    qwords: [u64; 16],
}

#[repr(C)]
#[derive(Copy, Clone)]
union BlockBuffer {
    bytes: [u8; BUFFER_SIZE * 128],
    blocks: [Block; BUFFER_SIZE],
}

const K: [u64; 80] = [
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
];

fn rotr(x: u64, n: u8) -> u64 {
    (x >> n) | (x << (64 - n))
}

fn ch(x: u64, y: u64, z: u64) -> u64 {
    (x & y) ^ (!x & z)
}

fn maj(x: u64, y: u64, z: u64) -> u64 {
    (x & y) ^ (x & z) ^ (y & z)
}

fn upper_sig0(x: u64) -> u64 {
    rotr(x, 28) ^ rotr(x, 34) ^ rotr(x, 39)
}

fn upper_sig1(x: u64) -> u64 {
    rotr(x, 14) ^ rotr(x, 18) ^ rotr(x, 41)
}

fn lower_sig0(x: u64) -> u64 {
    rotr(x, 1) ^ rotr(x, 8) ^ (x >> 7)
}

fn lower_sig1(x: u64) -> u64 {
    rotr(x, 19) ^ rotr(x, 61) ^ (x >> 6)
}

unsafe fn final_block(buffer: &mut BlockBuffer, bytes: usize, total_bits: u64) {
    let block: usize = (bytes - 1) / 128;
    let byte: usize = bytes % 128;

    if (byte < 112) & (byte > 0) {
        buffer.blocks[block].bytes[byte] = 0x80;
        for i in (byte + 1)..112 {
            buffer.blocks[block].bytes[i] = 0x00;
        }
        buffer.blocks[block].qwords[14] = 0;
        buffer.blocks[block].qwords[15] = total_bits.to_be();
    }

    else if (byte < 128) & (byte > 0) {
        buffer.blocks[block].bytes[byte] = 0x80;
        for i in (byte + 1)..128 {
            buffer.blocks[block].bytes[i] = 0x00;
        }
        if block < BUFFER_SIZE - 1 {
            for i in 0..14 {
                buffer.blocks[block + 1].qwords[i] = 0;
            }
            buffer.blocks[block + 1].qwords[14] = 0;
            buffer.blocks[block + 1].qwords[15] = total_bits.to_be();
        }
    }

    else if block < BUFFER_SIZE - 1 {
        for i in 0..14 {
            buffer.blocks[block + 1].qwords[i] = 0;
        }
        buffer.blocks[block + 1].bytes[0] = 0x80;
        buffer.blocks[block + 1].qwords[14] = 0;
        buffer.blocks[block + 1].qwords[15] = total_bits.to_be();
    }
}

unsafe fn add_padding(buffer: &mut BlockBuffer, total_bits: u64, add_eos: bool) {
    for i in 0..14 {
        buffer.blocks[0].qwords[i] = 0;
    }

    if add_eos {
        buffer.blocks[0].bytes[0] = 0x80;
    }

    buffer.blocks[0].qwords[14] = 0;
    buffer.blocks[0].qwords[15] = total_bits.to_be();
}

fn sha512_rounds(buffer: &mut BlockBuffer, state: &mut [u64; 8], iter: usize) {
    let (mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h, mut t1, mut t2): (u64, u64, u64, u64, u64, u64, u64, u64, u64, u64);
    let mut w: [u64; 80] = [0; 80];

    for block in 0..iter {
        for i in 0..16 {
            unsafe {
                w[i] = buffer.blocks[block].qwords[i].to_be();
            }
        }

        for i in 16..80 {
            w[i] = lower_sig1(w[i - 2]) + w[i - 7] + lower_sig0(w[i - 15]) + w[i - 16];
        }

        a = state[0];
        b = state[1];
        c = state[2];
        d = state[3];
        e = state[4];
        f = state[5];
        g = state[6];
        h = state[7];

        for i in 0..80 {
            t1 = h + upper_sig1(e) + ch(e, f, g) + K[i] + w[i];
            t2 = upper_sig0(a) + maj(a, b, c);
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

fn sha512(file_path: &str, iter: usize) -> std::io::Result<()> {
    let bytes: usize;
    let total_bits: u64;

    let mut f = File::open(file_path)?;

    let metadata = f.metadata()?;
    let total_bytes: usize = metadata.len() as usize;

    let mut buffers: BlockBuffer = BlockBuffer {
        bytes: [0; BUFFER_SIZE * 128],
    };

    let mut state: [u64; 8] = [
        0x6a09e667f3bcc908, 0xbb67ae8584caa73b, 
        0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1, 
        0x510e527fade682d1, 0x9b05688c2b3e6c1f, 
        0x1f83d9abfb41bd6b, 0x5be0cd19137e2179
    ];

    if total_bytes == 0 {
        unsafe {
            buffers.blocks[0].bytes[0] = 0x80;
        }
        for i in 1..128 {
            unsafe {
                buffers.blocks[0].bytes[i] = 0x00;
            }
        }
        sha512_rounds(&mut buffers, &mut state, 1);
        for i in 0..8 {
            print!("{:02x}", state[i]);
        }
        return Ok(());
    }

    total_bits = total_bytes as u64 * 8;

    for _ in 0..(total_bytes - 1) / (BUFFER_SIZE * 128 as usize) {
        unsafe {
            _ = f.read(&mut buffers.bytes[..])?;
        }
        for _ in 0..iter {
            sha512_rounds(&mut buffers, &mut state, BUFFER_SIZE);
        }
    }

    unsafe {
        bytes = f.read(&mut buffers.bytes[..])?;
    }
    unsafe {
        final_block(&mut buffers, bytes, total_bits);
    }

    if bytes > BUFFER_SIZE * 128 - 17 {
        sha512_rounds(&mut buffers, &mut state, BUFFER_SIZE);
        unsafe {
            add_padding(&mut buffers, total_bits, bytes % 128 == 0);
        }
        sha512_rounds(&mut buffers, &mut state, 1)
    }
    else {
        sha512_rounds(&mut buffers, &mut state, (bytes + 144) / 128);
    }

    println!();
    for i in 0..8 {
        print!("{:02x}", state[i]);
    }

    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 3 {
        println!("Usage: ./sha512_rs <File Path> <Iterations>");
        std::process::exit(0);
    }

    let file_path: &str = &args[1];

    let _ = sha512(file_path, args[2].parse::<usize>().unwrap());
}
