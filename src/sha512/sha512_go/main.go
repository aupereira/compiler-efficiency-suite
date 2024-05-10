package main

import (
	"encoding/binary"
	"fmt"
	"log"
	"os"
	"strconv"
)

const BufferSize uint64 = 1024

var K = [80]uint64{
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
}

func to_big_endian(n uint64) uint64 {
	return (n >> 56) | (n>>40)&0x000000000000FF00 |
		(n>>24)&0x0000000000FF0000 | (n>>8)&0x00000000FF000000 |
		(n<<8)&0x000000FF00000000 | (n<<24)&0x0000FF0000000000 |
		(n<<40)&0x00FF000000000000 | (n << 56)
}

func rotr(x uint64, n uint8) uint64 {
	return (x >> n) | (x << (64 - n))
}

func ch(x uint64, y uint64, z uint64) uint64 {
	return (x & y) ^ ((^x) & z)
}

func maj(x uint64, y uint64, z uint64) uint64 {
	return (x & y) ^ (x & z) ^ (y & z)
}

func upper_sig0(x uint64) uint64 {
	return rotr(x, 28) ^ rotr(x, 34) ^ rotr(x, 39)
}

func upper_sig1(x uint64) uint64 {
	return rotr(x, 14) ^ rotr(x, 18) ^ rotr(x, 41)
}

func lower_sig0(x uint64) uint64 {
	return rotr(x, 1) ^ rotr(x, 8) ^ (x >> 7)
}

func lower_sig1(x uint64) uint64 {
	return rotr(x, 19) ^ rotr(x, 61) ^ (x >> 6)
}

func final_block(buffer *[]uint64, bytes uint64, total_bits uint64) {
	block := (bytes - 1) / 128
	block_i := block * 16
	cbyte := bytes % 128
	cqword := cbyte / 8

	if cbyte < 112 && cbyte > 0 {
		(*buffer)[block_i+cqword] = (*buffer)[block_i+cqword] & ^(0xFF>>(8*(cbyte%8))) | 0x80<<(8*(cbyte%8))
		for i := cbyte/8 + 1; i < 14; i++ {
			(*buffer)[block*16+i] = 0
		}
		(*buffer)[block_i+14] = 0
		(*buffer)[block_i+15] = to_big_endian(total_bits)

	} else if cbyte < 128 && cbyte > 0 {
		(*buffer)[block_i+cqword] = (*buffer)[block_i+cqword] & ^(0xFF>>(8*(cbyte%8))) | 0x80<<(8*(cbyte%8))
		if cqword+1 == 15 {
			(*buffer)[block_i+15] = 0
		}

		if block < BufferSize-1 {
			block_i += 16
			for i := uint64(0); i < 14; i++ {
				(*buffer)[block_i+i] = 0
			}
			(*buffer)[block_i+14] = 0
			(*buffer)[block_i+15] = to_big_endian(total_bits)
		}

	} else if block < BufferSize-1 {
		block_i += 16
		(*buffer)[block_i] = 0x80
		for i := uint64(1); i < 14; i++ {
			(*buffer)[block_i+i] = 0
		}
		(*buffer)[block_i+14] = 0
		(*buffer)[block_i+15] = to_big_endian(total_bits)
	}
}

func add_padding(buffer *[]uint64, total_bits uint64, add_eos bool) {
	for i := 0; i < 14; i++ {
		(*buffer)[i] = 0
	}

	if add_eos {
		(*buffer)[0] = 0x80
	}

	(*buffer)[14] = 0
	(*buffer)[15] = to_big_endian(total_bits)
}

func sha512_rounds(buffer *[]uint64, state *[8]uint64, iter uint64) {
	var a, b, c, d, e, f, g, h, t1, t2 uint64
	var w [80]uint64 = [80]uint64{}

	var base uint64

	for block := uint64(0); block < iter; block++ {
		base = block * 16

		for i := uint64(0); i < 16; i++ {
			w[i] = to_big_endian((*buffer)[base+i])
		}

		for i := uint64(16); i < 80; i++ {
			w[i] = lower_sig1(w[i-2]) + w[i-7] + lower_sig0(w[i-15]) + w[i-16]
		}

		a = state[0]
		b = state[1]
		c = state[2]
		d = state[3]
		e = state[4]
		f = state[5]
		g = state[6]
		h = state[7]

		for i := 0; i < 80; i++ {
			t1 = h + upper_sig1(e) + ch(e, f, g) + K[i] + w[i]
			t2 = upper_sig0(a) + maj(a, b, c)
			h = g
			g = f
			f = e
			e = d + t1
			d = c
			c = b
			b = a
			a = t1 + t2
		}

		state[0] += a
		state[1] += b
		state[2] += c
		state[3] += d
		state[4] += e
		state[5] += f
		state[6] += g
		state[7] += h
	}
}

func sha512(file_path string, iter int) {
	var bytes int
	var r_error error

	f, err := os.Open(file_path)
	if err != nil {
		log.Fatal(err)
	}

	buffers := make([]byte, BufferSize*128)
	qbuffers := make([]uint64, BufferSize*16)

	state := [8]uint64{
		0x6a09e667f3bcc908, 0xbb67ae8584caa73b,
		0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1,
		0x510e527fade682d1, 0x9b05688c2b3e6c1f,
		0x1f83d9abfb41bd6b, 0x5be0cd19137e2179,
	}

	fi, err := f.Stat()
	if err != nil {
		log.Fatal(err)
	}
	total_bytes := uint64(fi.Size())
	total_bits := total_bytes * 8

	if total_bytes == 0 {
		qbuffers[0] = 0x80
		for i := 1; i < 16; i++ {
			buffers[i] = 0
		}
		sha512_rounds(&qbuffers, &state, 1)
		for i := 0; i < 8; i++ {
			fmt.Printf("%016x", state[i])
		}
		return
	}

	for _i := uint64(0); _i < (total_bytes-1)/(BufferSize*128); _i++ {
		bytes, r_error = f.Read(buffers)
		if r_error != nil {
			log.Fatal(r_error)
		}
		for _j := 0; _j < iter; _j++ {
			for i := 0; i < ((bytes-1)/8)+1; i++ {
				qbuffers[i] = binary.LittleEndian.Uint64(buffers[i*8 : (i*8)+8])
			}
			sha512_rounds(&qbuffers, &state, BufferSize)
		}
	}

	bytes, r_error = f.Read(buffers)
	if r_error != nil {
		log.Fatal(r_error)
	}
	for i := 0; i < ((bytes-1)/8)+1; i++ {
		qbuffers[i] = binary.LittleEndian.Uint64(buffers[i*8 : (i*8)+8])
	}
	final_block(&qbuffers, uint64(bytes), total_bits)

	if uint64(bytes) > BufferSize*128-17 {
		sha512_rounds(&qbuffers, &state, BufferSize)
		add_padding(&qbuffers, total_bits, bytes%128 == 0)
		sha512_rounds(&qbuffers, &state, 1)
	} else {
		sha512_rounds(&qbuffers, &state, (uint64(bytes)+144)/128)
	}

	for i := 0; i < 8; i++ {
		fmt.Printf("%016x", state[i])
	}

	defer f.Close()
}

func main() {
	if len(os.Args) != 3 {
		fmt.Println("Usage: ./sha512_go <File Path> <Iterations>")
		os.Exit(1)
	}

	arg2, err := strconv.ParseInt(os.Args[2], 10, 64)
	if err != nil {
		panic(err)
	}

	sha512(os.Args[1], int(arg2))
}
