#include <cstdlib>

void div90_main(int16_t poly[1440]);

int16_t poly[1440];
int main() {
  for (int i = 0; i < 1440; i++) {
    poly[i] = rand() % 4591 - 2295;
  }
  for (int t = 0; t < 1000000; t++) {
    poly[0] = rand() % 4591 - 2295;
    div90_main(poly);
  }
}
