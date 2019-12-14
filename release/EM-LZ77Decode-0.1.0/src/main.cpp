/**
 * @file    src/main.cpp
 * @section LICENCE
 *
 * Copyright (C) 2016-2019
 *   Djamal Belazzougui, Juha Karkkainen, Dominik Kempa, Simon J. Puglisi
 *
 * This file is part of EM-LZ77-Decode v0.1.0
 * Published in:
 *   Djamal Belazzougui, Juha Karkkainen,
 *   Dominik Kempa, Simon J. Puglisi:
 *   Lempel-Ziv Decoding in External Memory.
 *   Proc. 15th International Symposium on
 *   Experimental Algorithms (SEA), 2016.
 * 
 * Latest version: https://github.com/dkempa/em-lz77-decode
 * Contact: Dominik Kempa <dominik.kempa (at) gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 **/

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <cstring>
#include <getopt.h>
#include <unistd.h>
#include <omp.h>

#include "em_lz77decode_src/decode.h"


char *program_name;

void usage(int status) {
  printf(
"Usage: %s [OPTION]... FILE\n"
"Decode the LZ77 parsing stored in FILE.\n"
"\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -h, --help              display this help and exit\n"
"  -m, --mem=MEM           use MEM bytes of RAM for computation. Metric and IEC\n"
"                          suffixes are recognized, e.g., -l 10k, -l 1Mi, -l 3G\n"
"                          gives MEM = 10^4, 2^20, 3*10^6. Default: 3584Mi\n"
"  -o, --output=OUTFILE    specify the output file (default: FILE.decoded)\n",
    program_name);

  std::exit(status);
}

bool file_exists(std::string filename) {
  std::FILE *f = std::fopen(filename.c_str(), "r");
  bool ret = (f != NULL);
  if (f != NULL) std::fclose(f);

  return ret;
}

template<typename int_type>
bool parse_number(char *str, int_type *ret) {
  *ret = 0;
  std::uint64_t n_digits = 0;
  std::uint64_t str_len = std::strlen(str);
  while (n_digits < str_len && std::isdigit(str[n_digits])) {
    std::uint64_t digit = str[n_digits] - '0';
    *ret = (*ret) * 10 + digit;
    ++n_digits;
  }

  if (n_digits == 0)
    return false;

  std::uint64_t suffix_length = str_len - n_digits;
  if (suffix_length > 0) {
    if (suffix_length > 2)
      return false;

    for (std::uint64_t j = 0; j < suffix_length; ++j)
      str[n_digits + j] = std::tolower(str[n_digits + j]);
    if (suffix_length == 2 && str[n_digits + 1] != 'i')
      return false;

    switch(str[n_digits]) {
      case 'k':
        if (suffix_length == 1)
          *ret *= 1000;
        else
          *ret <<= 10;
        break;
      case 'm':
        if (suffix_length == 1)
          *ret *= 1000000;
        else
          *ret <<= 20;
        break;
      case 'g':
        if (suffix_length == 1)
          *ret *= 1000000000;
        else
          *ret <<= 30;
        break;
      case 't':
        if (suffix_length == 1)
          *ret *= 1000000000000;
        else
          *ret <<= 40;
        break;
      default:
        return false;
    }
  }

  return true;
}

int main(int argc, char **argv) {
  srand(time(0) + getpid());
  program_name = argv[0];

  static struct option long_options[] = {
    {"help",    no_argument,       NULL, 'h'},
    {"mem",     required_argument, NULL, 'm'},
    {"output",  required_argument, NULL, 'o'},
    {NULL, 0, NULL, 0}
  };

  std::uint64_t ram_use = 3584UL << 20;
  std::string out_fname("");

  // Parse command-line options.
  int c;
  while ((c = getopt_long(argc, argv, "hm:o:", long_options, NULL)) != -1) {
    switch(c) {
      case 'm':
        {
          bool ok = parse_number(optarg, &ram_use);
          if (!ok) {
            fprintf(stderr, "Error: parsing phrase length "
                "limit (%s) failed\n\n", optarg);
            usage(EXIT_FAILURE);
          }
          if (ram_use == 0) {
            fprintf(stderr, "Error: invalid RAM limit (%lu)\n\n", ram_use);
            usage(EXIT_FAILURE);
          }
          break;
        }
      case 'o':
        out_fname = std::string(optarg);
        break;
      case 'h':
        usage(EXIT_FAILURE);
      default:
        usage(EXIT_FAILURE);
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "Error: FILE not provided\n\n");
    usage(EXIT_FAILURE);
  }

  // Parse the text filename.
  std::string input_fname = std::string(argv[optind++]);
  if (optind < argc) {
    fprintf(stderr, "Warning: multiple input files provided. "
    "Only the first will be processed.\n");
  }

  // Set default output filename (if not provided).
  if (out_fname.empty())
    out_fname = input_fname + ".decoded";

  // Check if input exists.
  if (!file_exists(input_fname)) {
    fprintf(stderr, "Error: input file (%s) does not exist\n\n",
        input_fname.c_str());
    usage(EXIT_FAILURE);
  }

  if (file_exists(out_fname)) {

    // Output file exists, should we proceed?
    char *line = NULL;
    std::uint64_t buflen = 0;
    std::int64_t len = 0L;

    do {
      printf("Output file (%s) exists. Overwrite? [y/n]: ",
          out_fname.c_str());
      if ((len = getline(&line, &buflen, stdin)) == -1) {
        printf("\nError: failed to read answer\n\n");
        std::fflush(stdout);
        usage(EXIT_FAILURE);
      }
    } while (len != 2 || (line[0] != 'y' && line[0] != 'n'));

    if (line[0] == 'n') {
      free(line);
      std::exit(EXIT_FAILURE);
    }
    free(line);
  }

  // Run the algorithm decoding LZ77.
  lz77decode(input_fname, out_fname, ram_use);
}
