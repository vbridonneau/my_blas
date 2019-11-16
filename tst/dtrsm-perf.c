#include "dtrsm.h"
#include "algonum.h"
#include <stdio.h>

void test_result() {

}

int main(int argc, char **argv) {
	char *str[] = {"OK", "NOK"};

	// int r1 = !!testone_dtrsm( my_dtrsm, 1000, 1000, 1 );
	int r1 = !!testall_dtrsm( my_dtrsm );
	printf("%s\n", str[r1]);

	test_result();
	return 0;
}
