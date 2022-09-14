#include <stdio.h>
#include <time.h>
long long cpucycles(void) {
  unsigned long long result;
  __asm__ volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax" : "=a" (result) ::  "%rdx");
  return result;
}

int main(){
	long long c1, c2;
	time_t clk1, clk2;
	int j,i;

	c1 = cpucycles();
	j=0;
	for(i=0; i<100000000; i++){
	}
	c2 = cpucycles();
	printf("%Ld %Ld %Ld\n",c1, c2, c2-c1);


	clk1 = clock();
	j=0;
	for(i=0; i<100000000; i++){
	}
	clk2 = clock();
	printf("%ld %ld %ld\n",clk1, clk2, clk2-clk1);

	return 0;
}
