/**
 * @file    src/em_lz77decode_src/decode.hpp
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

#ifndef __EM_LZ77DECODE_SRC_DECODE_HPP_INCLUDED
#define __EM_LZ77DECODE_SRC_DECODE_HPP_INCLUDED

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>

#include "io/async_multi_stream_writer.hpp"
#include "io/async_stream_reader.hpp"
#include "io/async_stream_writer.hpp"
#include "io/async_vbyte_stream_reader.hpp"


namespace em_lz77decode_private {

template<typename value_type>
inline void append_vbyte_encoding(
    std::uint8_t *buf,
    std::uint64_t &filled,
    value_type value) {
  std::uint64_t ptr = filled;
  while (value > 127) {
    buf[ptr++] = ((value & 0x7f) | 0x80);
    value >>= 7;
  }
  buf[ptr++] = value;
  filled = ptr;
}

inline std::uint64_t vbyte_encoding_size(std::uint64_t value) {
  std::uint64_t ret = 1UL;
  while (value > 127) { ++ret; value >>= 7; }
  return ret;
}


template<typename reader_type, typename return_type>
inline return_type decode_vbyte_from_stream(reader_type *reader) {
  return_type result = 0;
  std::uint64_t offset = 0;
  std::uint8_t next_byte = reader->read();
  while (next_byte & 0x80) {
    result |= (((return_type)next_byte & 0x7F) << offset);
    offset += 7;
    next_byte = reader->read();
  }
  result |= ((return_type)next_byte << offset);
  return result;
}

//===========================================================================
// For any phrase (pos, len) starting at position 'beg' write the triple
// (pos, beg, len) into file corresponding to the segment containing position
// 'pos' (i.e., phrase source).
//===========================================================================
template<std::uint64_t alignment_unit_log>
bool permute_triples_by_source(
    std::string parsing_filename,
    std::string output_filename,
    std::uint64_t max_segment_size,
    std::uint64_t disk_space_use,
    std::uint64_t in_buf_ram,
    std::uint64_t &text_length,
    std::uint64_t &n_phrases,
    std::uint64_t &parsing_file_ptr,
    std::uint64_t &total_io_volume) {

  static const std::uint64_t alignment_unit = (1UL << alignment_unit_log);
  static const std::uint64_t alignment_unit_mask = alignment_unit - 1;

  // Set initial values.
  std::uint64_t text_len = text_length;
  std::uint64_t n_phr = n_phrases;
  std::uint64_t n_segments = 0;
  std::uint64_t segments_reach = 0;
  std::uint64_t dbg = 0;

  // Specify reader/write type.
  typedef async_vbyte_stream_reader<std::int64_t> reader_type;
  typedef async_multi_stream_writer<std::uint8_t> multifile_writer_type;

  // Create reader/writer.
  reader_type *parsing_reader =
    new reader_type(parsing_filename, parsing_file_ptr, in_buf_ram);
  multifile_writer_type *writer = new multifile_writer_type();

  // Allocate small auxiliary buffer.
  static const std::uint64_t bufsize = (1UL << 10);
  std::uint8_t *buf = new std::uint8_t[bufsize];

  std::vector<std::uint64_t> recent_phrase_end;

  std::uint64_t parsing_file_size =
    utils::file_size(parsing_filename);

  fprintf(stderr, "  Total parsing data decoded so far: %.2Lf%%\n",
      (100.L * parsing_file_ptr) / parsing_file_size);
  fprintf(stderr, "  Total text decoded so far: %lu (%.2LfMiB)\n",
      text_length, 1.L * text_length / (1L << 20));

  // Scan the parsing and place phrases in buckets.
  fprintf(stderr, "\r  Permute phrases by source: ");
  long double permute_start = utils::wclock();

  // The upper bound on the peak disk space usage (in
  // addition to parsing) when processing current part.
  std::uint64_t cur_disk_space_use_upper_bound = text_len;

  // chunks_encoding_size[j][i] tells the size of
  // chunks encoding sampled from segment i with the
  // destination in segment j. We deliberatelly invert
  // the matrix to reduce the number of cache misses.
  std::vector<std::vector<std::uint64_t> > chunks_encoding_size;
  std::vector<std::uint64_t> triple_encoding_size;  // for phrases in the source in given segment

  // Used as an upper bound on the vbyte encoding
  // size of 'position' component of the chunk header.
  std::uint64_t vbyte_encoding_size_of_segment_pos =
    vbyte_encoding_size(max_segment_size);

  while (!parsing_reader->empty()) {

    // Print progress message.
    ++dbg;
    if (dbg == (1L << 20)) {
      dbg = 0;
      long double elapsed = utils::wclock() - permute_start;
      std::uint64_t io_vol = parsing_reader->bytes_read() + writer->bytes_written();
      fprintf(stderr, "\r  Permute phrases by source: time = %.2Lfs, I/O = %.2LfMiB/s",
          elapsed, (1.L * io_vol / (1 << 20)) / elapsed);
    }

    if (cur_disk_space_use_upper_bound >= disk_space_use) {

      // Compute a much better upper bound on the peak disk usage.
      // First, compute suffix sum of triple_encoding_size.
      std::vector<std::uint64_t>
        sum_of_triple_encoding_size_suffix_sum(n_segments);
      for (std::uint64_t segment_id = n_segments;
          segment_id > 0; --segment_id) {
        sum_of_triple_encoding_size_suffix_sum[segment_id - 1] =
          triple_encoding_size[segment_id - 1];
        if (segment_id != n_segments)
          sum_of_triple_encoding_size_suffix_sum[segment_id - 1] +=
            sum_of_triple_encoding_size_suffix_sum[segment_id];
      }

      // Go through segments and for each compute much
      // better upper bound on the peak disk usage.
      std::uint64_t better_upper_bound_on_disk_usage = 0UL;
      std::uint64_t tmp = 0;  // reduces complexity O(n_segments^3) -> O(n_segments^2).

      for (std::uint64_t segment_id = 0; segment_id < n_segments; ++segment_id) {
        std::uint64_t this_segment_beg = segment_id * max_segment_size;
        std::uint64_t this_segment_end =
          std::min(this_segment_beg + max_segment_size, text_len);
        std::uint64_t this_segment_size = this_segment_end - this_segment_beg;

        // This is the space occupied at all times
        // when processing the current segment.
        std::uint64_t this_segment_fixed_space =
          this_segment_beg;  // text decoded so far
        if (segment_id + 1 != n_segments)  // triples encoding for segments to the right
          this_segment_fixed_space +=
            sum_of_triple_encoding_size_suffix_sum[segment_id + 1];

        if (segment_id == 0) {

          // For first segment we compute it in O(n_segments^2).
          tmp = 0;
          for (std::uint64_t t = segment_id + 1; t < n_segments; ++t)
            for (std::uint64_t tt = 0; tt < segment_id; ++tt)
              tmp += chunks_encoding_size[t][tt];
        } else {

          // For other segment we modify the value from
          // previous iteration in O(n_segments) time.
          for (std::uint64_t t = segment_id + 1; t < n_segments; ++t)
            tmp += chunks_encoding_size[t][segment_id - 1];
          for (std::uint64_t tt = 0; tt < segment_id - 1; ++tt)
            tmp -= chunks_encoding_size[segment_id][tt];
        }

        this_segment_fixed_space += tmp;

        // Space for chunks coming *into* current segment.
        std::uint64_t this_segment_chunks_to = 0;
        for (std::uint64_t t = 0; t < segment_id; ++t)
          this_segment_chunks_to += chunks_encoding_size[segment_id][t];

        // Space for chunks coming *from* current segment.
        std::uint64_t this_segment_chunks_from = 0UL;
        for (std::uint64_t t = segment_id + 1; t < n_segments; ++t)
          this_segment_chunks_from += chunks_encoding_size[t][segment_id];

        // Get the space for the triples file of the current segment.
        std::uint64_t this_segment_triples_file_size =
          triple_encoding_size[segment_id];

        // We are now ready to compute the actual peak disk space usage
        // during processing of the current segment, taking into account
        // in which order which file is deleted / expanded.
        // Note: it would be possibly to further optimize this estimation
        // if we kept the triples corresponding to segment in two files,
        // one containing the triples with the position in the segment,
        // and the other with triples in other segmnets (to the right).
        // This would help here, but increase the number of buffers
        // necessary in the phrase permutation step.
        std::uint64_t this_segment_peak_disk_usage = 0UL;
        if (segment_id + 1 == n_segments) {
          this_segment_peak_disk_usage =
            this_segment_fixed_space +
            std::max(this_segment_triples_file_size +
                this_segment_chunks_to, this_segment_size);
        } else {
          this_segment_peak_disk_usage =
            this_segment_fixed_space +
            this_segment_triples_file_size +
            std::max(this_segment_chunks_from +
                this_segment_size, this_segment_chunks_to);
        }

        // Compare the result for the current peak.
        better_upper_bound_on_disk_usage =
          std::max(better_upper_bound_on_disk_usage,
              this_segment_peak_disk_usage);
      }

      if (better_upper_bound_on_disk_usage >= disk_space_use) break;
      else cur_disk_space_use_upper_bound =
        better_upper_bound_on_disk_usage;
    }

    // Extract next phrase.
    std::int64_t pos = parsing_reader->read();                  // source of current phrase
    std::uint64_t len = (std::uint64_t)parsing_reader->read();  // length of current phrase
    std::uint64_t beg = text_len;                               // beginning of current phrase

    // Update text length and parsing size.
    text_len += std::max(1UL, len);
    ++n_phr;

    // Create new writer(s) if the new phrase extends
    // the region covered by existing segments.
    while (segments_reach < text_len) {
      std::string new_seg_filename = output_filename +
        ".seg" + utils::intToStr(n_segments++) + ".triples";
      writer->add_file(new_seg_filename);
      recent_phrase_end.push_back(0);
      for (std::uint64_t j = 0; j + 1 < n_segments; ++j)
        chunks_encoding_size[j].push_back(0);
      std::vector<std::uint64_t> vec(n_segments, 0);
      chunks_encoding_size.push_back(vec);
      triple_encoding_size.push_back(0);
      segments_reach += max_segment_size;
    }

    // Add triple (pos, beg, len) to the file corresponding
    // to segment in which 'pos' lies (if len > 0). All values
    // are vbyte-encoded. Furthermore, rather than storing beg
    // we store the distance to the end of the most recently
    // added phrase (in segment 'segment_id'). Also, rather
    // than storing the pos values, we store the vbyte encoding
    // of pos - dest_segment_beg, where dest_segment_beg is the
    // beginning of the destination segment.
    if (!len) {
      std::uint64_t dest_segment_id = beg / max_segment_size;
 
      std::uint64_t filled = 0;
      append_vbyte_encoding(buf,
          filled, beg - recent_phrase_end[dest_segment_id]);
      append_vbyte_encoding(buf, filled, pos);
      append_vbyte_encoding(buf, filled, len);
      writer->write_to_ith_file(dest_segment_id, buf, filled);
      recent_phrase_end[dest_segment_id] = beg + 1;

      triple_encoding_size[dest_segment_id] += filled;
      cur_disk_space_use_upper_bound += filled;
    } else {
      while (true) {
        std::uint64_t l_segment_left =
          alignment_unit - (pos & alignment_unit_mask);
        std::uint64_t r_segment_left =
          alignment_unit - (beg & alignment_unit_mask);
        std::uint64_t curlen =
          std::min(std::min(l_segment_left, r_segment_left), len);
        std::uint64_t beg_segment_id = beg / max_segment_size;
        std::uint64_t dest_segment_id = pos / max_segment_size;
        std::uint64_t dest_segment_beg = dest_segment_id * max_segment_size;

        std::uint64_t filled = 0;
        append_vbyte_encoding(buf, filled,
            beg - recent_phrase_end[dest_segment_id]);
        append_vbyte_encoding(buf, filled,
            pos - dest_segment_beg);

        // Used to get the size of vbyte encoding of curlen.
        std::uint64_t filled_snapshot = filled;
        append_vbyte_encoding(buf, filled, curlen);
        writer->write_to_ith_file(dest_segment_id, buf, filled);

        triple_encoding_size[dest_segment_id] += filled;
        chunks_encoding_size[beg_segment_id][dest_segment_id] +=
          curlen +
          vbyte_encoding_size_of_segment_pos +
          (filled - filled_snapshot);
        cur_disk_space_use_upper_bound +=
          filled +
          curlen +
          vbyte_encoding_size_of_segment_pos +
          (filled - filled_snapshot);

        recent_phrase_end[dest_segment_id] = beg + curlen;
        if (curlen == len) break;
        else {
          len -= curlen;
          pos += curlen;
          beg += curlen;
        }
      }
    }
  }

  bool end_of_parsing_file = (parsing_reader->empty());

  // Update I/O volume.
  std::uint64_t io_vol =
    parsing_reader->bytes_read() +
    writer->bytes_written();
  total_io_volume += io_vol;

  // Print summary.
  long double permute_time = utils::wclock() - permute_start;
  fprintf(stderr, "\r  Permute phrases by source: 100.0%%, "
      "time = %.2Lfs, I/O = %.2LfMiB/s\n",
      permute_time, ((1.L * io_vol) / (1 << 20)) / permute_time);
  fprintf(stderr, "  Processing this part will decode "
      "%.1LfMiB of text (%.2Lf%% of parsing)\n",
      1.L * (text_len - text_length) / (1L << 20),
      (100.L * parsing_reader->bytes_read()) / parsing_file_size);
  fprintf(stderr, "  Last part = %s\n",
      end_of_parsing_file ? "TRUE" : "FALSE");

  // Update output arguments of the function.
  parsing_file_ptr += parsing_reader->bytes_read();
  text_length = text_len;
  n_phrases = n_phr;

  // Clean up.
  delete parsing_reader;
  delete writer;
  delete[] buf;

  return end_of_parsing_file;
}

//==============================================================================
// Place the text chunks from file (if it exists)
// in the correct place of the current segment.
//==============================================================================
void load_segment_chunks(std::uint8_t *segment, std::string chunk_filename,
    std::uint64_t in_buf_ram, std::uint64_t &total_io_volume) {
  typedef async_stream_reader<std::uint8_t> chunk_reader_type;
  chunk_reader_type *chunk_reader =
    new chunk_reader_type(chunk_filename,
        in_buf_ram, in_buf_ram / (2UL << 20));

  std::uint64_t dbg = 0;
  std::uint64_t chunk_file_bytes = utils::file_size(chunk_filename);

  fprintf(stderr, "\r    Load segment chunks: ");
  long double load_chunks_start = utils::wclock();
  while (!chunk_reader->empty()) {
     // Print progress message.
    if (dbg >= (4L << 20)) {
      long double elapsed = utils::wclock() - load_chunks_start;
      std::uint64_t io_vol = chunk_reader->bytes_read();
      fprintf(stderr, "\r    Load segment chunks: %.1Lf%%, "
          "time = %.2Lfs, I/O = %.2LfMiB/s",
          (1.L * io_vol * 100.L) / chunk_file_bytes,
          elapsed, (1.L * io_vol / (1L << 20)) / elapsed);
      dbg = 0;
    }

    // Deserialize a single chunk from stream.
    std::uint64_t chunk_dest =
      decode_vbyte_from_stream<chunk_reader_type, std::uint64_t>(chunk_reader);
    std::uint64_t chunk_len =
      decode_vbyte_from_stream<chunk_reader_type, std::uint64_t>(chunk_reader);
    for (std::uint64_t t = 0; t < chunk_len; ++t)
      segment[chunk_dest + t] = chunk_reader->read();

    // Update number of read bytes and I/O volume.
    dbg += chunk_len;
  }

  // Print summary.
  long double load_chunks_time =
    utils::wclock() - load_chunks_start;
  fprintf(stderr, "\r    Load segment chunks: "
      "100.0%%, time = %.2Lfs, I/O = %.2LfMiB/s\n",
      load_chunks_time,
      (1.L * chunk_file_bytes / (1L << 20)) / load_chunks_time);

  // Clean up.
  total_io_volume += chunk_reader->bytes_read();
  delete chunk_reader;
}

//=============================================================================
// Process all triples with the source in the current segment. Processing
// consits of (1) updating the current segment, (2) sampling the symbols
// from the current segment and distributing into buckets containing chunks.
//=============================================================================
void process_triples(
    std::uint8_t *segment,
    std::uint64_t segment_id,
    std::uint64_t text_length,
    std::uint64_t max_segment_size,
    std::string output_filename,
    std::uint64_t in_buf_ram,
    std::uint64_t out_buf_ram,
    std::uint64_t &total_io_volume,
    std::uint64_t precomputed_prefix_length,
    std::uint64_t uncomputed_suffix_length) {

  typedef async_vbyte_stream_reader<std::uint64_t> triple_reader_type;
  std::uint64_t segment_beg = segment_id * max_segment_size;
  std::uint64_t segment_end =
    std::min(segment_beg + max_segment_size, text_length);
  std::uint64_t next_segment_end =
    std::min(text_length, segment_end + max_segment_size);

  std::string triples_filename =
    output_filename + ".seg" +
    utils::intToStr(segment_id) + ".triples";
  triple_reader_type *triple_reader =
    new triple_reader_type(triples_filename, 0UL, in_buf_ram);

  std::uint64_t previous_phrase_end = 0;
  bool triple_avail = false;
  std::uint64_t absolute_beg = 0;
  std::uint64_t relative_pos = 0;
  std::uint64_t len = 0;
  std::uint64_t dbg = 0;
  std::uint64_t io_vol = 0;

  if (!triple_reader->empty()) {
    absolute_beg = previous_phrase_end + triple_reader->read();
    relative_pos = triple_reader->read();
    len = triple_reader->read();
    previous_phrase_end = absolute_beg + std::max(1UL, len);
    triple_avail = true;
  }

  if (triple_avail) {
    fprintf(stderr, "\r    Selfdecode segment: ");
    long double selfdecode_start = utils::wclock();

    while (triple_avail) {

      // Print progress message.
      ++dbg;
      if (dbg == (1L << 20)) {
        long double elapsed = utils::wclock() - selfdecode_start;
        io_vol = triple_reader->bytes_read();
        fprintf(stderr, "\r    Selfdecode segment: time = %.2Lfs, "
            "I/O = %.2LfMiB/s",
            elapsed, (1.L * io_vol / (1L << 20)) / elapsed);
        dbg = 0;
      }

      if (absolute_beg >= segment_end) break;

      std::uint64_t relative_beg = absolute_beg - segment_beg;
      if (!len) segment[relative_beg] = (std::uint8_t)relative_pos;
      else {
        for (std::uint64_t t = 0; t < len; ++t)
          segment[relative_beg + t] = segment[relative_pos + t];
      }

      // Read another triple.
      if (triple_reader->empty()) triple_avail = false;
      else {
        absolute_beg = previous_phrase_end + triple_reader->read();
        relative_pos = triple_reader->read();
        len = triple_reader->read();
        previous_phrase_end = absolute_beg + std::max(1UL, len);
      }
    }

    // Update I/O volume.
    io_vol = triple_reader->bytes_read();
    total_io_volume += io_vol;

    // Print summary.
    long double selfdecode_time = utils::wclock() - selfdecode_start;
    fprintf(stderr, "\r    Selfdecode segment: "
        "time = %.2Lfs, I/O = %.2LfMiB/s\n",
        selfdecode_time,
        ((1.L * io_vol) / (1L << 20)) / selfdecode_time);
  }

  // Allocate auxiliary buffer.
  static const std::uint64_t bufsize = (1UL << 10);
  std::uint8_t *buf = new std::uint8_t[bufsize];

  // Store the state of the computation.
  std::uint64_t previous_phrase_end_copy = previous_phrase_end;
  bool triple_avail_copy = triple_avail;
  std::uint64_t absolute_beg_copy = absolute_beg;
  std::uint64_t relative_pos_copy = relative_pos;
  std::uint64_t len_copy = len;
  std::uint64_t triple_file_pos = triple_reader->bytes_read();

  // Distribtue chunks with the source in current segment.
  // Optimization: only distribute the ones with
  // absolute_beg - absolute_pos >= max_segment_size.
  if (triple_avail) {
    fprintf(stderr, "\r    Distribute chunks: ");
    long double distribute_start = utils::wclock();

    typedef async_stream_writer<std::uint8_t> writer_type;
    writer_type *chunk_writer = NULL;
    std::uint64_t prev_dest_segment_id = (1UL << 60);

    dbg = 0;
    std::uint64_t write_io_volume = 0;
    std::uint64_t triple_reader_bytes_subtract =
      triple_reader->bytes_read();

    while (triple_avail) {
      ++dbg;
      if (dbg == (1L << 20)) {
        dbg = 0;
        long double elapsed = utils::wclock() - distribute_start;
        io_vol =
          (triple_reader->bytes_read() - triple_reader_bytes_subtract) +
          write_io_volume;
        fprintf(stderr, "\r    Distribute chunks: "
            "time = %.2Lfs, I/O = %.2LfMiB/s",
            elapsed, (1.L * io_vol / (1L << 20)) / elapsed);
      }

      std::uint64_t absolute_pos = segment_beg + relative_pos;
      if (absolute_beg - absolute_pos >= max_segment_size) {
        // Open the new writer if necessary.
        std::uint64_t dest_segment_id = absolute_beg / max_segment_size;
        if (dest_segment_id != prev_dest_segment_id) {
          prev_dest_segment_id = dest_segment_id;
          std::string filename =
            output_filename + ".seg" +
            utils::intToStr(dest_segment_id) + ".chunks";

          if (chunk_writer != NULL) delete chunk_writer;
          chunk_writer = new writer_type(filename, out_buf_ram,
              out_buf_ram / (2UL << 20), "a");
        }

        // Process next triple.
        std::uint64_t dest_segment_beg = dest_segment_id * max_segment_size;
        std::uint64_t relative_beg = absolute_beg - dest_segment_beg;
        dest_segment_id -= (segment_id + 1);

        // Serialize the chunk into the chunk file
        // corresponding to segment dest_segment_id.
        std::uint64_t filled = 0;
        append_vbyte_encoding(buf, filled, relative_beg);
        append_vbyte_encoding(buf, filled, len);

        std::uint64_t tocopy = std::min(bufsize - filled, len);
        std::copy(segment + relative_pos,
            segment + relative_pos + tocopy, buf + filled);
        filled += tocopy;
        chunk_writer->write(buf, filled);
        write_io_volume += filled;

        if (len > tocopy) {
          len -= tocopy;
          relative_pos += tocopy;
          chunk_writer->write(segment + relative_pos, len);
          write_io_volume += len;
        }
      }

      // Fetch next triple.
      if (triple_reader->empty()) triple_avail = false;
      else {
        absolute_beg = previous_phrase_end + triple_reader->read();
        relative_pos = triple_reader->read();
        len = triple_reader->read();
        previous_phrase_end = absolute_beg + std::max(1UL, len);
      }
    }

    // Update I/O volume.
    io_vol =
      (triple_reader->bytes_read() - triple_reader_bytes_subtract) +
      write_io_volume;
    total_io_volume += io_vol;

    // Print summary.
    long double distribute_time = utils::wclock() - distribute_start;
    fprintf(stderr, "\r    Distribute chunks: "
        "time = %.2Lfs, I/O = %.2LfMiB/s\n",
        distribute_time,
        (1.L * io_vol / (1L << 20)) / distribute_time);

    // Clean up.
    if (chunk_writer != NULL)
      delete chunk_writer;
  }

  // Append uncomputed segment suffix to output file.
  if (uncomputed_suffix_length > 0) {

    // Print initial message.
    fprintf(stderr, "    Write segment to disk: ");
    long double write_start = utils::wclock();
    std::FILE *f = utils::file_open(output_filename, "a");
    utils::write_to_file(segment + precomputed_prefix_length,
        uncomputed_suffix_length, f);

    // Clean up.
    std::fclose(f);

    // Update I/O volume.
    total_io_volume += uncomputed_suffix_length;

    // Print summary.
    long double write_time = utils::wclock() - write_start;
    fprintf(stderr, "%.2Lfs, I/O = %.2LfMiB/s\n",
        write_time,
        (1.L * uncomputed_suffix_length / (1L << 20)) / write_time);
  }

  delete[] buf;
  delete triple_reader;
  triple_reader = NULL;

  fprintf(stderr, "    Decode the nearby phrases for next segment: ");
  long double nearby_decode_start = utils::wclock();

  // Scan one more time the section of the triples
  // file that contains phrsaes with
  // absolute_beg - absolute_pos < max_segment_size
  // and decode then using window symbols.
  triple_reader =
    new triple_reader_type(triples_filename, triple_file_pos);
  triple_avail = triple_avail_copy;
  previous_phrase_end = previous_phrase_end_copy;
  relative_pos = relative_pos_copy;
  len = len_copy;
  absolute_beg = absolute_beg_copy;

  while (triple_avail) {

    // Process next triple.
    if (absolute_beg >= next_segment_end) break;
    std::uint64_t absolute_pos = segment_beg + relative_pos;
    if (absolute_beg - absolute_pos < max_segment_size) {

      // Decode symbols [relative_beg..relative_beg+len) of the next segment
      // from symbols [relative_pos..relative_pos+len) in the current segment.
      std::uint64_t relative_beg = absolute_beg - segment_end;
      for (std::uint64_t j = 0; j < len; ++j)
        segment[relative_beg + j] = segment[relative_pos + j];
    }

    // Fetch next triple.
    if (triple_reader->empty()) triple_avail = false;
    else {
      absolute_beg = previous_phrase_end + triple_reader->read();
      relative_pos = triple_reader->read();
      len = triple_reader->read();
      previous_phrase_end = absolute_beg + std::max(1UL, len);
    }
  }

  // Print summary.
  long double nearby_decode_time =
    utils::wclock() - nearby_decode_start;
  fprintf(stderr, "%.2Lfs\n", nearby_decode_time);

  // Clean up.
  delete triple_reader;
}

void decode(
    std::string parsing_filename,
    std::string output_filename,
    std::uint64_t ram_use_incl_buffers,
    std::uint64_t disk_space_use) {

  // Turn paths absolute.
  parsing_filename = utils::absolute_path(parsing_filename);
  output_filename = utils::absolute_path(output_filename);
  utils::drop_disk_pages(parsing_filename);

  // Unit of alignment, used to simplify handling long phrases.
  static const std::uint64_t alignment_unit_log = 22L;
  static const std::uint64_t in_buf_ram = (16UL << 20);
  static const std::uint64_t out_buf_ram = (16UL << 20);
  std::uint64_t ram_use = ram_use_incl_buffers - in_buf_ram - out_buf_ram;
  std::uint64_t parsing_size_bytes = utils::file_size(parsing_filename);
  std::uint64_t max_segment_size_units =
    std::max(1UL, ram_use >> alignment_unit_log);
  std::uint64_t max_segment_size =
    (max_segment_size_units << alignment_unit_log);
  std::uint64_t text_length = 0;
  std::uint64_t n_phrases = 0;
  std::uint64_t parsing_file_ptr = 0;
  std::uint64_t n_parts = 0;
  std::uint64_t total_io_volume = 0;

  // Print initial info.
  fprintf(stderr, "Input filename = %s\n", parsing_filename.c_str());
  fprintf(stderr, "Output filename = %s\n", output_filename.c_str());
  fprintf(stderr, "Input size = %ld (%.1LfMiB)\n",
      parsing_size_bytes, 1.L * parsing_size_bytes / (1L << 20));
  fprintf(stderr, "Max disk use = %ld (%.1LfMiB)\n",
      disk_space_use, 1.L * disk_space_use / (1L << 20));
  fprintf(stderr, "RAM use = %ld (%.1LfMiB)\n",
      ram_use_incl_buffers, 1.L * ram_use_incl_buffers / (1L << 20));
  fprintf(stderr, "RAM use (excl. buffer) = %ld (%.1LfMiB)\n",
      ram_use, 1.L * ram_use / (1L << 20));
  fprintf(stderr, "Max segment size = %ld (%.1LfMiB)\n",
      max_segment_size, (1.L * max_segment_size) / (1L << 20));

  // Delete the output file, if exists.
  if (utils::file_exists(output_filename))
     utils::file_delete(output_filename);

  long double start = utils::wclock();

  bool end_of_parsing = false;
  while (!end_of_parsing) {
    ++n_parts;
    fprintf(stderr, "\nProcess part %ld:\n", n_parts);

    // Permute triples by source (note:
    // permute_triples_by_source updates text_length).
    std::uint64_t old_text_length = text_length;
    end_of_parsing =
      permute_triples_by_source<alignment_unit_log>(parsing_filename,
        output_filename, max_segment_size, disk_space_use, in_buf_ram,
        text_length, n_phrases, parsing_file_ptr, total_io_volume);

    // Process segments left to right.
    std::uint8_t *cur_segment = new std::uint8_t[max_segment_size];
    std::uint64_t n_segments =
      (text_length + max_segment_size - 1) / max_segment_size;

    for (std::uint64_t segment_id = 0; segment_id < n_segments; ++segment_id) {
      std::uint64_t segment_beg = segment_id * max_segment_size;
      std::uint64_t segment_end =
        std::min(segment_beg + max_segment_size, text_length);
      std::uint64_t segment_size = segment_end - segment_beg;

      std::uint64_t precomputed_prefix_length = 0UL;
      if (old_text_length > segment_beg)
        precomputed_prefix_length =
          std::min(old_text_length - segment_beg, segment_size);
      std::uint64_t uncomputed_suffix_length =
        segment_size - precomputed_prefix_length;

      std::string chunk_filename =
        output_filename + ".seg" + utils::intToStr(segment_id) + ".chunks";
      std::string triples_filename = output_filename + ".seg" +
        utils::intToStr(segment_id) + ".triples";
      bool nonempty_chunks_file =
        (utils::file_exists(chunk_filename) &&
         utils::file_size(chunk_filename) > 0);
      bool nonempty_triples_file =
        (utils::file_exists(triples_filename) &&
         utils::file_size(triples_filename) > 0);

      if (nonempty_chunks_file || nonempty_triples_file) {
        fprintf(stderr, "  Decode segment %ld/%ld [%ld..%ld):\n",
            segment_id + 1, n_segments, segment_beg, segment_end);

        // Load precomputed segment prefix from file.
        if (precomputed_prefix_length > 0) {
          fprintf(stderr, "    Read segment from disk: ");
          long double read_start = utils::wclock();
          std::uint64_t toread =
            std::min(old_text_length - segment_beg, segment_size);
          utils::read_at_offset(cur_segment,
              segment_beg, toread, output_filename);

          // Update I/O volume.
          total_io_volume += toread;

          // Print summary.
          long double read_time = utils::wclock() - read_start;
          fprintf(stderr, "%.2Lfs, I/O = %.2LfMiB/s\n",
              read_time, (1.L * toread / (1L << 20)) / read_time);
        }

        // Load chunks inside current segment.
        if (nonempty_chunks_file) {
          load_segment_chunks(cur_segment,
              chunk_filename, in_buf_ram, total_io_volume);
          utils::file_delete(chunk_filename);
        }

        // Process phrases with the source in the current segment.
        if (nonempty_triples_file) {
          process_triples(cur_segment, segment_id, text_length,
              max_segment_size, output_filename, in_buf_ram,
              out_buf_ram, total_io_volume, precomputed_prefix_length,
              uncomputed_suffix_length);
          utils::file_delete(triples_filename);
        }
      }

      if (utils::file_exists(chunk_filename))
        utils::file_delete(chunk_filename);
      if (utils::file_exists(triples_filename))
        utils::file_delete(triples_filename);
    }

    // Clean up.
    delete[] cur_segment;
  }

  // Print summary.
  long double total_time = utils::wclock() - start;
  fprintf(stderr, "\n\nComputation finished. Summary:\n");
  fprintf(stderr, "  I/O volume = %lu (%.2LfB/B)\n",
      total_io_volume, 1.L * total_io_volume / text_length);
  fprintf(stderr, "  number of decoded phrases: %ld\n", n_phrases);
  fprintf(stderr, "  average phrase length = %.2Lf\n",
      (1.L * text_length) / n_phrases);
  fprintf(stderr, "  length of decoded text: %ld (%.2LfMiB)\n",
      text_length, (1.L * text_length) / (1L << 20));
  fprintf(stderr, "  elapsed time: %.2Lfs (%.4Lfs/MiB)\n",
      total_time, total_time / ((1.L * text_length) / (1L << 20)));
  fprintf(stderr, "  speed: %.2LfMiB/s\n",
      ((1.L * text_length) / (1L << 20)) / total_time);
}

}  // namespace em_lz77decode_private

void lz77decode(
    std::string input_filename,
    std::string out_filename,
    std::uint64_t ram_use,
    std::uint64_t disk_space_use = (1UL << 60)) {
  em_lz77decode_private::decode(input_filename,
      out_filename, ram_use, disk_space_use);
}

#endif  // __EM_LZ77DECODE_SRC_DECODE_HPP_INCLUDED
