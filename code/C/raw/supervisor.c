#include <stdio.h>

int main(){
  system("git add docs");
  system("git commit -m \"website update\" ");
  system("git push");
  return 0;
}
