#include "print.h"

#ifdef FOUND_PYBIND11
#include <pybind11/embed.h>

namespace py = pybind11;

void printInTqdm(const char* str)
{
  py::module_ tqdm = py::module_::import("tqdm");
  py::object tqdm_write = tqdm.attr("tqdm").attr("write");
  tqdm_write(str);
}
#else
void printInTqdm(const char* str)
{
  std::cout << str;
}
#endif

namespace std {
  std::ostream tqdm_cout(new callback_streambuf(printInTqdm));
}

void printInTqdmFormat(const char* format, ...)
{
  va_list args;
  va_start(args, format);
  char buffer[1024];
  vsprintf(buffer, format, args);
  va_end(args);
  printInTqdm(buffer);
}
