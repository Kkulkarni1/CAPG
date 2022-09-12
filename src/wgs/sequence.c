#include "sequence.h"

size_t _data_size = 8 * sizeof(data_t);

//extern data_t *sequence_alloc(size_t nchar, sequence_opt *so);
extern data_t read_char(sequence *s, sequence_opt *so, size_t i);
//extern void write_char(sequence *s, sequence_opt *so, size_t i, data_t in);
//extern size_t seqlen(sequence *s, sequence_opt *so);
//extern size_t seqbytes(sequence *s, sequence_opt *so);
