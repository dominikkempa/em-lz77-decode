EM-LZ77-Decode - external memory decoder of LZ77 parsing
========================================================


Description
-----------

This package contains implementation of the external memory algorithm
for decoding of LZ77 parsing. The algorithm is described in the paper

  Djamal Belazzougui, Juha Karkkainen, Dominik Kempa, Simon J. Puglisi:
  Lempel-Ziv Decoding in External Memory.
  15th International Symposium on Experimental Algorithms (SEA), 2016.

The latest version of the algorithm is available from
https://github.com/dkempa/em-lz77-decode.



Requirements
------------

EM-LZ77-Decode has no external dependencies (libraries, cmake, etc).
It only requires:
- g++ compiler supporting the -std=c++0x flag (all modern versions)
- A 64-bit operating system. The current version has been tested
  on Linux/PC.



Compilation and usage
---------------------

The package contains a single Makefile in the main directory. Type
'make' to build the executable. For usage instructions, run the
program without any arguments.

Example
~~~~~~~

The simplest usage of EM-LZ77-Decode is as follows. Suppose the
vbyte-encoded (see below) LZ77 parsing is located in
/data/input.txt.lz77.vbyte. Then, to decode the text type:

  $ ./decode_lz77 /data/input.txt.lz77.vbyte

This will write the text to /data/input.txt.lz77.vbyte.decoded.
By default, the algorithm uses 3.5GiB of RAM for computation. A
more advanced usage is demonstrated below.

  $ ./decode_lz77 ./input.txt.lz77.vbyte -m 8gi -o ../input.txt

Explanation:
- The 'm' flag allows specifying the amount of RAM used during the
  computation (in bytes). In this example, the RAM limit is set to
  8gi = 8 * 2^30 bytes.
- The 'o' flag allows specifying the location and filename of the
  output text. The default location and filename is the same
  as input text, with the appended ".decoded" suffix.

Notes:
- The argument of the 'm' flag (RAM used during the computation)
  can be specified either explicitly or using common suffixes such
  as K, M, G, T, Ki, Mi, Gi, Ti which respectively correspond to
  multipliers: 10^3, 10^6, 10^9, 10^12, 2^10, 2^20 2^30, 2^40. Suffix
  names are not case-sensitive, e.g., Ti = ti, k = K.
- The flags specifying RAM usage, output filename, etc. can be
  given in any order.
- Filenames passed as command-line arguments can be given as absolute,
  relative, and common (such as $HOME) paths, e.g., ../input.txt and
  ~/data/input.txt are valid paths, see also example above.
- EM-LZ77-Decode expects the pairs in the input LZ77 parsing to be
  stored using vbyte encoding. A standard encoding of pairs using
  fixed-width integers (produced by, e.g., LZ77 parsing algorithm
  EM-LPF) can be converted into the vbyte version using the tool
  located in the ./tools directory.



Troubleshooting
---------------

1. I am getting an error about the exceeded number of opened files.

Solution: The error is caused by the operating system imposing a limit
on the maximum number of files opened by a program. The limit can be
increased with the "ulimit -n newlimit" command. However, in Linux the
limit cannot be increased beyond the so-called "hard limit", which is
usually only few times larger. Furthermore, this is a temporary
solution that needs to repeated every time a new session is
started. To increase the limits permanently, edit (as a root) the file
/etc/security/limits.conf and add the following lines at the end
(including the asterisks):

* soft nofile 128000
* hard nofile 128000

This increases the limit to 128000 (use larger values if necessary).
The new limits apply (check with ulimit -n) after starting new
session.



Limitations
-----------

- Current implementation supports only decoding of texts over byte
  alphabet.
- At present the only limitation in the usage of the algorithm is the
  need to ensure that the limit for the number of opened files in the
  system is sufficiently large to prevent the above error. This
  technical shortcoming will be eliminated in the future versions.



Terms of use
------------

EM-LZ77-Decode is released under the MIT/X11 license. See the file
LICENCE for more details.

If you use this code, please cite the paper mentioned above and
publish the URL from which you downloaded the code.
