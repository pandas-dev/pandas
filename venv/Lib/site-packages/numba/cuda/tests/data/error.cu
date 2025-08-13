extern "C" __device__
int bar(int* out, int a) {
  // Explicitly placed to generate an error
  SYNTAX ERROR
  *out = a * 2;
  return 0;
}
