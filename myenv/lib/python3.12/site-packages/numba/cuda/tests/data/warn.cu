extern "C" __device__
int bar(int* out, int a) {
  // Explicitly placed to generate a warning for testing the NVRTC program log
  int unused;
  *out = a * 2;
  return 0;
}
