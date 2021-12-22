#include <stdlib.h>

#include "nuc.h"

/**
 * XY used penultimate 2 bits to encode 4 nucleotides.
 * IUPAC re-encodes using 4 bits to represent 17 characters, ignoring U as not
 *   distinct from T.
 *
 * char ASCII   binary xy_t iupac_t    binary
 *    A    65  1000001    0       1  00000001
 *    C    67  1000011    1       2  00000010
 *    G    71  1000111    3       4  00000100
 *    T    84  1010100    2       8  00001000
 *    U    85  1010101            8  00001000
 *    R    82  1010010            5  00000101
 *    Y    89  1011001           10  00001010
 *    S    83  1010011            6  00000110
 *    W    87  1010111            9  00001001
 *    K    75  1001011           12  00001100
 *    M    77  1001101            3  00000011
 *    B    66  1000010           14  00001110
 *    D    68  1000100           13  00001101
 *    H    72  1001000           11  00001011
 *    V    86  1010110            7  00000111
 *    N    78  1001110    -      15  00001111
 *    X    88  1011000            0  00000000
 */

/**
 * Number of nucleotides consistent with each IUPAC symbol. [NOT USED]
 */
const unsigned char popcnt[] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};

/**
 * Convert xy_t to char for human consumption.
 */
unsigned char const xy_to_char[NUM_NUCLEOTIDES] = {'A', 'C', 'T', 'G'};

/**
 * Convert std_t to char for human consumption.
 */
unsigned char const std_to_char[NUM_NUCLEOTIDES] = {'A', 'C', 'G', 'T'};

/**
 * Convert iupac_t to char for human consumption.
 */
unsigned char const iupac_to_char[NUM_IUPAC_SYMBOLS] = {
	'-', 'A', 'C', 'M', 'G', 'R', 'S',
	'V', 'T', 'W', 'Y', 'H', 'K',
	'D', 'B', 'N'
};

/**
 * Convert xy_t to iupac_t.
 */
iupac_t const xy_to_iupac[NUM_NUCLEOTIDES] = {
	IUPAC_A,
	IUPAC_C,
	IUPAC_T,
	IUPAC_G
};

/**
 * Convert xy_t to reverse complement xy_t.
 */
xy_t const xy_to_rc[NUM_NUCLEOTIDES] = {
	XY_T,
	XY_G,
	XY_A,
	XY_C
};

/**
 * Convert iupac_t to reverse complement iupac_t.
 */
iupac_t const iupac_to_rc[NUM_IUPAC_SYMBOLS] = {
	15, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 0
};

/**
 * Convert iupac_t to xy_t: only meaningful for A, C, G, T.
 */
xy_t const iupac_to_xy[NUM_IUPAC_SYMBOLS] = {
	0, XY_A, XY_C, 0, XY_G, 0, 0, 0, XY_T, 0, 0, 0, 0, 0, 0, 0
};

/**
 * Convert iupac_t to standard nucleotide order A=0, C=1, G=2, T=3, which is
 * NOT xy_t.
 */
unsigned char const iupac_to_std[NUM_IUPAC_SYMBOLS] = {
	0, STD_A, STD_C, 0, STD_G, 0, 0, 0, STD_T, 0, 0, 0, 0, 0, 0, 0
};

/**
 * Standard order: A, C, G, T
 * XY order: A, C, T, G
 */
unsigned char const xy_to_std[NUM_NUCLEOTIDES] = {0, 1, 3, 2};
xy_t const std_to_xy[NUM_NUCLEOTIDES] = {0, 1, 3, 2};

/**
 * Convert char-encoded nucleotide - MIN_NUCLEOTIDE_ASCII to iupac_t
 * type.  Handles 'A' thru 'Y' and converts unknown to 0 aka 'X'.
 */
iupac_t const nuc_to_iupac[NUCLEOTIDE_ALPHABET_SIZE] = {
	1,14,2,13,0,0,4,11,0,0,12,0,3,15,0,0,0,5,6,8,8,7,9,0,10
	/*                      1                    2
        0  1 2  3 4 5 6  7 8 9  0 1 2  3 4 5 6 7 8 9 0 1 2 3  4 */
};

/**
 * Convert char-encoded nucleotide - MIN_NUCLEOTIDE_ASCII to xy_t
 * type.  Handles 'A' thr 'Y' but converts unknown to 0 aka 'A'.
 */
xy_t const nuc_to_xy[NUCLEOTIDE_ALPHABET_SIZE] = {
	0,0,1,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,0
	/*                  1                   2
        0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 */
};

/*
const uint8_t _xy_bits = 2;
const uint8_t _xy_mask = 3;
const uint8_t _iupac_bits = 4;
const uint8_t _iupac_mask = 15;
*/
#define _xy_bits	2U
#define _xy_mask	3U
#define _iupac_bits	4U
#define _iupac_mask	15U

sequence_opt _xy_sequence_opt = {_xy_bits, _xy_mask};
sequence_opt _iupac_sequence_opt = {_iupac_bits, _iupac_mask};

/**
 * Validate human-readable nucleotide characters as IUPAC symbols or one
 * of standard nucleotides: A, C, G, T.
 */
extern int valid_iupac(char *c);
extern int valid_nucleotide(char *c);
extern xy_t char_to_xy(char c);
extern iupac_t char_to_iupac(char c);
extern data_t char_to_data(char c, int encoding);

extern unsigned int number_nucleotide(iupac_t c);
extern const iupac_t *nucleotide_list(iupac_t c);

extern sequence_opt *nuc_sequence_opt(int encoding);
extern int fwrite_nuc_sequence(FILE *fp, sequence *s, int encoding);
extern int fwrite_nuc_segment(FILE *fp, sequence *s, int encoding, size_t sidx, size_t eidx);
extern void fwrite_xy_sequence(FILE *fp, sequence *s, size_t sidx, size_t eidx);
extern void fwrite_iupac_sequence(FILE *fp, sequence *s, size_t sidx, size_t eidx);
extern data_t get_nuc(sequence *s, int encoding, size_t i);

