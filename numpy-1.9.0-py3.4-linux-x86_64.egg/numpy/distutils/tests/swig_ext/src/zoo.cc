#include "zoo.h"
#include <cstdio>
#include <cstring>

Zoo::Zoo()
{
  n = 0;
}

void Zoo::shut_up(char *animal)
{
  if (n < 10) {
    strcpy(animals[n], animal);
    n++;
  }
}

void Zoo::display()
{
  int i;
  for(i = 0; i < n; i++)
    printf("%s\n", animals[i]);
}
