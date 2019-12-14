#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <getopt.h>
#include <unistd.h>

#include "stream.h"
#include "utils.h"
#include "uint40.h"


template<typename writer_type>
inline void encode_and_write_vbyte(writer_type *writer, long value) {
  while (value > 127) {
    writer->write((value & 0x7f) | 0x80);
    value >>= 7;
  }
  writer->write(value);
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s FILE\n"
        "Convert parsing in FILE to vbyte encoding (pos and "
        "len) and write to FILE.vbyte.\n", argv[0]);
    std::exit(EXIT_FAILURE);
  }

  std::string parsing_filename = argv[1];
  std::string output_filename = parsing_filename + ".vbyte";

  long n_phrases = utils::file_size(parsing_filename) / 10;

  typedef stream_reader<uint40> reader_type;
  typedef stream_writer<unsigned char> writer_type;

  reader_type *reader = new reader_type(parsing_filename);
  writer_type *writer = new writer_type(output_filename);

  for (long i = 0; i < n_phrases; ++i) {
    long pos = reader->read();
    long len = reader->read();
    encode_and_write_vbyte(writer, pos);
    encode_and_write_vbyte(writer, len);
  }

  delete reader;
  delete writer;
}
