import numpy
import random
import sys
 
def random_message(length):
    msg = []
    for i in range(length):
        letter = random.choice([0,1])
        msg.append(letter)
    return numpy.array(msg, dtype = int)
 
def encode(m, g):
    enc = numpy.dot(m, g)%2
    return enc
 
def get_syndrome(m, h):
    dec = numpy.dot(h, m)%2
    return dec
 
def noise(m, error, bit):
    noisy = numpy.copy(m)
    cont = 0
    for i in range(len(noisy)):
        e = random.random()
        if e <= error:
            if noisy[i] == 0:
                noisy[i] = 1
                cont +=1
            else:
                noisy[i] = 0
                cont +=1
            if cont == bit:
                break
    return noisy
 
def get_error(m):
    n = ''
    for i in range(len(m)):
        n += str(m[i])
        #print "n: ", n
    #n = n[::-1]
    #print "n: ", n

    n = (int(n,2)) -1
    #print "n: ", n

    return n
 
def correct(noisy, n):
    if noisy[n] == 0:
        noisy[n] = 1
    else:# noisy[n] == 1:
        noisy[n] = 0
    #print noisy
    return noisy
 
def hamming(length, n, error, bit):
    g =  numpy.array([[1, 0, 0, 0, 0, 1, 1],[0, 1, 0, 0, 1, 0, 1],[0, 0, 1, 0, 1, 1, 0],[0, 0, 0, 1, 1, 1, 1]])
    h = numpy.array([[0, 0, 0, 1, 1, 1, 1],[0, 1, 1, 0, 0, 1, 1],[1, 0, 1, 0, 1, 0, 1],])
    msg = random_message(length)
    print "original message: ", msg
    msg_chunks = []
    while len(msg) >= 4:
        msg_chunk = msg[0:4]
        corrected = numpy.copy(msg_chunk)
        for i in range(n):
            
            enc = encode(msg_chunk, g)
            print "original message chunk: ", msg_chunk
            print "encoded msg: ", enc
            noisy = noise(enc, error, bit)
            print 'noisy mesage: ', noisy
            syndrome = get_syndrome(noisy, h)
            print 'syndrome: ', syndrome
            error_bit = get_error(syndrome)
            if error_bit >= 0:
                corrected = correct(noisy, error_bit)
        print 'corrected: ', corrected
        msg_chunks.append(corrected[0:4])
        msg = msg[4:]
    decoded = []
    for msg_chunk in msg_chunks:
        decoded = decoded + list(msg_chunk)
    print 'decoded: ', decoded
            
        
if __name__ == '__main__':
    length = int(raw_input('Length: '))
    n = int(raw_input('Repetitions: '))
    error = float(raw_input('Error percentage: '))
    bit = int(raw_input('Bits with error: '))
    hamming(length,n,error, bit)