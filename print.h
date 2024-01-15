#ifndef PRINT_H
#define PRINT_H

#include <iostream>
#include <stdarg.h>
#include <functional>

#define printf(...) printInTqdmFormat(__VA_ARGS__)
#define cout tqdm_cout

void printInTqdm(const char* str);
void printInTqdmFormat(const char* format, ...);

class callback_streambuf : public std::streambuf {
public:
  callback_streambuf(std::function<void(char const*)> callback) : callback(callback) {}

protected:
  std::streamsize xsputn(char_type const* s, std::streamsize count) {
    std::string str(s, s + count);
    size_t pos = str.find('\n');
    if (pos != std::string::npos) {
      std::string line = str.substr(0, pos + 1);
      callback(line.c_str());
      str = str.substr(pos + 1);
    }
    return count;
  }

  int sync() {
    if (!str.empty()) {
      callback(str.c_str());
      str.clear();
    }
    return 0;
  }

private:
  std::function<void(char const*)> callback;
  std::string str;
};

namespace std {
  std::ostream extern tqdm_cout;
}

#endif