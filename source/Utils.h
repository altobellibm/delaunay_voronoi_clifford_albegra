#ifndef UTILS_H
#define UTILS_H

namespace Utils{
	
	// https://en.wikipedia.org/wiki/Hamming_weight
	// http://stackoverflow.com/questions/109023/how-to-count-the-number-of-set-bits-in-a-32-bit-integer
	int HammingWeight(unsigned int i)
	{
		i = i - ((i >> 1) & 0x55555555);
		i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
		return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
	}

};

#endif 
