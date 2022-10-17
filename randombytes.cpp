#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include "randombytes.hpp"
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
//#include <sys/syscall.h>

void randombytes(uint8_t *out, size_t outlen) {
  static int fd = -1;
  ssize_t ret;

  while(fd == -1) {
    fd = open("/dev/urandom", O_RDONLY);
    if(fd == -1 && errno == EINTR)
      continue;
    else if(fd == -1)
      abort();
  }

  while(outlen > 0) {
    ret = read(fd, out, outlen);
    if(ret == -1 && errno == EINTR)
      continue;
    else if(ret == -1)
      abort();

    out += ret;
    outlen -= ret;
  }
}
