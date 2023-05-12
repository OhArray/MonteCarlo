import numpy as np
import random
import math             


class BitString:
    """
    Simple class to implement a string of bits
    """
    def __init__(self, string):
        self.string = string

        self.N = len(self)
        self.n_dim = 2**self.N
        
    def __str__(self):
        bitstring = "".join(map(str, self.string))
        return str(bitstring)
        
    def __len__(self):
        return len(self.string)
    
    def __getitem__(self,i):
        return int(self.string[i])

    '''
    def __eq__(self, other):
        if len(self.string) != len(other.string):
            return False
        for i in range(len(self.string)):
            if self.string != other.string:
                return False
        return True
    '''
    
    def flip(self, bit):
        state = self.string[bit]
        if state == 0:
            self.string[bit] = 1
        else:
            self.string[bit] = 0
    
    def set_string(self, bitstring):
        self.string = bitstring
        
    def on(self):
        return self.string.count(1)
    
    def off(self):
        return self.string.count(0)
    
    def int(self):
        sum = 0
        upper = len(self)
        for i in range(upper):
            curr = upper - i - 1
            sum += self.string[i] * (2**(curr))
        return sum
    
    def set_int(self, decimal, digits=None):
        binary = []
        count = 0
        while decimal > 0:
            bit = decimal % 2
            binary.append(bit)
            decimal = (decimal - bit) // 2
            count += 1
        if digits == None:
            digits = self.N
        if digits != None:
            for i in range(count, digits):
                binary.append(0)
        binary.reverse()
        self.string = binary

    def set_magnetization(self, M = 0):
        randomlist = random.sample(range(0, self.N), M)
        for i in randomlist:
            self.string[i] = 1

    def get_magnetization(self, M = 0):
        return np.sum(2*self.array()-1)

    def array(self):
        x = np.array(self.string)
        y = x.astype(np.int)
        return y
