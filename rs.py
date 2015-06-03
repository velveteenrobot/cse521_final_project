import numpy
import random
import sys
import datetime



class RSError(Exception):
    pass

gf_exp = [1] * 512
gf_log = [0] * 256
x = 1



for i in range(1, 255):
    x <<= 1
    if x & 0x100:
        x ^= 0x11d
    gf_exp[i] = x
    gf_log[x] = i
for i in range(255, 512):
    gf_exp[i] = gf_exp[i - 255]

def gf_mul(x, y):
    if x == 0 or y == 0:
        return 0
    return gf_exp[gf_log[x] + gf_log[y]]

def gf_div(x, y):
    if y == 0:
        raise ZeroDivisionError()
    if x == 0:
        return 0
    return gf_exp[gf_log[x] + 255 - gf_log[y]]

def gf_poly_scale(p, x):
    return [gf_mul(p[i], x) for i in range(0, len(p))]

def gf_poly_add(p, q):
    r = [0] * max(len(p), len(q))
    for i in range(0, len(p)):
        r[i + len(r) - len(p)] = p[i]
    for i in range(0, len(q)):
        r[i + len(r) - len(q)] ^= q[i]
    return r

def gf_poly_mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for j in range(0, len(q)):
        for i in range(0, len(p)):
            r[i + j] ^= gf_mul(p[i], q[j])
    return r

def gf_poly_eval(p, x):
    y = p[0]
    for i in range(1, len(p)):
        y = gf_mul(y, x) ^ p[i]
    return y

def rs_generator_poly(nsym):
    g = [1]
    for i in range(0, nsym):
        g = gf_poly_mul(g, [1, gf_exp[i]])
    return g

def rs_encode_msg(msg_in, nsym):
    if len(msg_in) + nsym > 255:
        raise ValueError("message too long")
    gen = rs_generator_poly(nsym)
    msg_out = msg_in + [0] * nsym
    for i in range(0, len(msg_in)):
        coef = msg_out[i]
        if coef != 0:
            for j in range(0, len(gen)):
                msg_out[i + j] ^= gf_mul(gen[j], coef)
    msg_out[:len(msg_in)] = msg_in
    return msg_out

def rs_calc_syndromes(msg, nsym):
    synd = [0] * nsym
    for i in range(0, nsym):
        synd[i] = gf_poly_eval(msg, gf_exp[i])
    return synd

def rs_correct_errata(msg, synd, pos):
    # calculate error locator polynomial
    q = [1]
    for i in range(0, len(pos)):
        x = gf_exp[len(msg) - 1 - pos[i]]
        q = gf_poly_mul(q, [x, 1])
    # calculate error evaluator polynomial
    p = synd[0:len(pos)]
    p.reverse()
    p = gf_poly_mul(p, q)
    p = p[len(p) - len(pos):len(p)]
    # formal derivative of error locator eliminates even terms
    q = q[len(q) & 1:len(q):2]
    # compute corrections
    for i in range(0, len(pos)):
        x = gf_exp[pos[i] + 256 - len(msg)]
        y = gf_poly_eval(p, x)
        z = gf_poly_eval(q, gf_mul(x, x))
        msg[pos[i]] ^= gf_div(y, gf_mul(x, z))


def rs_find_errors(synd, nmess):
    # find error locator polynomial with Berlekamp-Massey algorithm
    err_poly = [1]
    old_poly = [1]
    for i in range(0, len(synd)):
        old_poly.append(0)
        delta = synd[i]
        for j in range(1, len(err_poly)):
            delta ^= gf_mul(err_poly[len(err_poly) - 1 - j], synd[i - j])
        if delta != 0:
            if len(old_poly) > len(err_poly):
                new_poly = gf_poly_scale(old_poly, delta)
                old_poly = gf_poly_scale(err_poly, gf_div(1, delta))
                err_poly = new_poly
            err_poly = gf_poly_add(err_poly, gf_poly_scale(old_poly, delta))
    errs = len(err_poly) - 1
    if errs * 2 > len(synd):
        raise RSError("Too many errors to correct")
    # find zeros of error polynomial
    err_pos = []
    for i in range(0, nmess):
        if gf_poly_eval(err_poly, gf_exp[255 - i]) == 0:
            err_pos.append(nmess - 1 - i)
    if len(err_pos) != errs:
        return None    # couldn't find error locations
    return err_pos

def rs_forney_syndromes(synd, pos, nmess):
    fsynd = list(synd)      # make a copy
    for i in range(0, len(pos)):
        x = gf_exp[nmess - 1 - pos[i]]
        for i in range(0, len(fsynd) - 1):
            fsynd[i] = gf_mul(fsynd[i], x) ^ fsynd[i + 1]
        fsynd.pop()
    return fsynd

def rs_correct_msg(msg_in, nsym):
    if len(msg_in) > 255:
        raise ValueError("message too long")
    msg_out = list(msg_in)     # copy of message
    # find erasures
    erase_pos = []
    for i in range(0, len(msg_out)):
        if msg_out[i] < 0:
            msg_out[i] = 0
            erase_pos.append(i)
    if len(erase_pos) > nsym:
        raise RSError("Too many erasures to correct")
    synd = rs_calc_syndromes(msg_out, nsym)
    if max(synd) == 0:
        return msg_out  # no errors
    fsynd = rs_forney_syndromes(synd, erase_pos, len(msg_out))
    err_pos = rs_find_errors(fsynd, len(msg_out))
    if err_pos is None:
        raise RSError("Could not locate error")
    rs_correct_errata(msg_out, synd, erase_pos + err_pos)
    synd = rs_calc_syndromes(msg_out, nsym)
    if max(synd) > 0:
        raise RSError("Could not correct message")
    return msg_out[:-nsym]



class RS(object):
    total_time = 0
    def __init__(self, nsym=10):
        self.nsym = nsym

    def encode(self, data):
        start = datetime.datetime.now()
        if isinstance(data, str):
            data = [ord(ch) for ch in data]
        if not isinstance(data, list):
            data = list(data)
        chunk_size = 255 - self.nsym
        enc = []
        for i in range(0, len(data), chunk_size):
            chunk = data[i:i+chunk_size]
            enc.extend(rs_encode_msg(chunk, self.nsym))
        end = datetime.datetime.now()
        self.total_time += (end - start).total_seconds()
        return enc
    
    def decode(self, data):
        start = datetime.datetime.now()
        if isinstance(data, str):
            data = [ord(ch) for ch in data]
        dec = []
        for i in range(0, len(data), 255):
            chunk = data[i:i+255]
            dec2 = rs_correct_msg(chunk, self.nsym)
            if len(dec2) == len(chunk):
                # no errors, discard ECC
                dec2 = dec2[:-self.nsym]
            dec.extend(dec2)
        end = datetime.datetime.now()
        self.total_time += (end - start).total_seconds()
        return dec

    def add_noise(self, msg, error, bit):
        noisy = numpy.copy(msg)
        cont = 0
        for i in range(len(noisy)):
            e = random.random()
            if e <= error:
                noisy[i] += 1
                cont+= 1
                if cont == bit:
                    break
        return noisy

    def binary_chunk_to_dec(self, data, chunk_size):
        dec_chunk_list = []
        while len(data) >= chunk_size:
            chunk = data[:chunk_size]
            chunk = ''.join(map(str,chunk))
            chunk_dec = int(chunk, 2)
            dec_chunk_list.append(chunk_dec)
            data = data[chunk_size:]
        return dec_chunk_list

    def dec_to_binary_chunk(self, data, chunk_size):
        result = []
        for num in data:
            binary_num = bin(num)[2:]
            if len(binary_num) < chunk_size:
                diff = chunk_size - len(binary_num)
                binary_num = ''.join(['0']*diff) + binary_num
            binary_array = []
            for char in binary_num:
                binary_array.append(int(char))

            result = result + binary_array 
        return result

if __name__ == '__main__':
    length = int(raw_input('Length: '))
    n = int(raw_input('Repetitions: '))
    error = float(raw_input('Error percentage: '))
    bit = int(raw_input('Bits with error: '))
    parity_bits = int(raw_input('Parity bits: '))
    block_length = int(raw_input('Block length: '))
    
    msg = [1] * length
    rs = RS(parity_bits)
    a = rs.binary_chunk_to_dec(msg, block_length)
    for i in range(n):
        b = rs.encode(a)
        c = rs.add_noise(b, error, bit)
        d = rs.decode(c)

    print "Total time: " + str(rs.total_time/n)


